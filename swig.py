import ismrmrd
import logging
import traceback
import constants
import numpy as np
import ctypes
from utils import svd_compress, atleast_nd, resize
from ismrmrdtools import coils
from mrdhelper import update_img_header_from_raw

import sys
sys.path.append('/opt/code/bart/python')
from bart import bart

def groups(iterable, predicate):
    group = []
    for item in iterable:
        group.append(item)

        if predicate(item):
            yield group
            group = []

def conditionalGroups(iterable, predicateAccept, predicateFinish):
    group = []
    try:
        for item in iterable:
            if item is None:
                break

            if predicateAccept(item):
                group.append(item)

            if predicateFinish(item):
                yield group
                group = []
    except Exception as e:
        logging.error(traceback.format_exc())
        iterable.send_logging(constants.MRD_LOGGING_ERROR, traceback.format_exc())
    finally:
        iterable.send_close()

def process(connection, config, metadata):
    logging.info("Config: \n%s", config)
    logging.info("Metadata: \n%s", metadata)

    image_index_offset = 0

    # Discard phase correction lines and accumulate lines until "ACQ_LAST_IN_MEASUREMENT" is set
    for group in conditionalGroups(connection, lambda acq: not acq.is_flag_set(ismrmrd.ACQ_IS_PHASECORR_DATA), lambda acq: acq.is_flag_set(ismrmrd.ACQ_LAST_IN_MEASUREMENT)):
        images = process_group(group, config, metadata, image_index_offset)
        image_index_offset += len(images)
        connection.send_image(images)

def process_group(group, config, metadata, image_index_offset):

    data = [acquisition.data for acquisition in group]
    data = np.stack(data, axis=-1)

    user_int = [acquisition.user_int for acquisition in group]
    user_int = np.stack(user_int, axis=-1)

    physiology_time_stamp = [acquisition.physiology_time_stamp for acquisition in group]
    physiology_time_stamp = np.stack(physiology_time_stamp, axis=-1)

    acquisition_time_stamp = [acquisition.acquisition_time_stamp for acquisition in group]
    acquisition_time_stamp = np.stack(acquisition_time_stamp, axis=-1)

    sets = [acquisition.idx.set for acquisition in group]
    number_of_vencs = max(sets) + 1
    echo_spacing = getattr(metadata.sequenceParameters, "echo_spacing[0]", [0])[0]

    # Reformat data from [cha col lin] -> [col lin cha]
    data = np.moveaxis(data, 0, -1)

    # Perform coil compression if number of coils exceedes 8
    if data.shape[-1] > 8:
        data = svd_compress(data)

    # Reformat data from [col lin cha] -> [1 col lin cha] to match BART convention
    data = np.expand_dims(data, 0)
    _, ncol, nlin, ncha = data.shape
    logging.info("Raw data is size %s" % (data.shape,))

    # Extract spoke angles from header and convert to float
    angles = user_int[6, :].astype(np.float64)*np.pi/(2**14)

    traj = np.zeros((3,ncol,nlin))
    for venc in range(number_of_vencs):
        spokes_in_venc = [i for i, x in enumerate(sets) if x == venc]
        spokes_in_venc = np.array(spokes_in_venc, dtype=int)
        angles_in_venc = angles[spokes_in_venc]

        # Use 75 maximally spaced spokes for gradient delay estimation
        sorted_idx = np.argsort(angles_in_venc)
        delay_est_idx = np.linspace(0, len(sorted_idx) - 1, 75, dtype=int)

        # Calculate trajectory from selected angles only
        sorted_angles = angles_in_venc[sorted_idx[delay_est_idx]]
        sorted_traj = bart(1, f'traj -x{ncol} -r -c', C=sorted_angles)

        # Perform gradient delay estimation
        # Rosenzweig et. al. Simple Auto-Calibrated Gradient Delay Estimation From Few Spokes Using Radial Intersections, MRM 2018
        delay = bart(1, 'estdelay -R', sorted_traj, data[:, :, spokes_in_venc[sorted_idx[delay_est_idx]], :])
        delay_str = ":".join(map(str, delay.real))

        # Calculate gradient delay corrected trajectory for all angles
        traj_venc = bart(1, f'traj -x{ncol} -r -q{delay_str}, -c -O', C=angles_in_venc)
        traj[:, :, spokes_in_venc] = traj_venc.real

    # Reconstruct a temporal averaged image for coil sensitivity estimation
    coil_images = bart(1, f'nufft -d{ncol}:{ncol}:1 -i -t', traj, data)

    # BART output is [x y z cha] -> ismrmrd-python-tools input is [cha y x]
    coil_images = np.squeeze(coil_images).transpose(2, 1, 0)

    # Calculate coil sensitivity maps using the iterative Inati method from Gadgetron;
    # Inati SJ, Hansen MS, Kellman P. A solution to the phase problem in adaptive coil combination.
    sens,_ = coils.calculate_csm_inati_iter(coil_images)

    # Revert back to BART convention [cha y x] -> [x y z cha]
    sens = np.expand_dims(sens.transpose(2, 1, 0), 2)

    # Identify cardiac triggers
    ecg_time_stamp = physiology_time_stamp[0, :].astype(np.int16)
    first_in_beat = np.concatenate(([0], np.flatnonzero(np.diff(ecg_time_stamp)<0)))+1
    spokes_per_beat = np.diff(first_in_beat)
    number_of_beats = spokes_per_beat.size
    logging.info(f'Spokes per beat: {spokes_per_beat}')

    # Calculate heartbeat statistics
    # Time stamps are in ticks of 2.5 ms
    RR = np.diff(acquisition_time_stamp[first_in_beat])*2.5
    RR_avg = int(np.round(np.mean(RR)))
    RR_std = int(np.round(np.std(RR)))
    logging.info(f'RR {RR_avg} +/- {RR_std}; {number_of_beats} heartbeats')

    window_width = 3
    number_of_phases = int(np.median(spokes_per_beat)/window_width)

    logging.info(f'Window width: {window_width}')
    logging.info(f'Number of phases: {number_of_phases}')

    # Calculate trigger times (in ticks) to be written back into the image header
    triggerTime = np.linspace(0, (RR_avg-(window_width/2)*echo_spacing)/2.5, number_of_phases*2, dtype=int)

    # Perform binning as described paper
    C = np.arange(number_of_phases)/(number_of_phases+1)
    cardiac_phase = np.empty([window_width*number_of_beats, number_of_vencs, number_of_phases], dtype=int)
    for venc in range(number_of_vencs):
        spokes_in_venc = [i for i, x in enumerate(sets) if x == venc]
        for beat in range(number_of_beats):
            spokes_in_beat = np.intersect1d(np.arange(first_in_beat[beat], first_in_beat[beat+1], dtype=int), spokes_in_venc)
            number_of_spokes_in_beat = spokes_in_beat.size
            cyclic_spoke_axis = np.tile(np.arange(number_of_spokes_in_beat, dtype=int), 3)

            spoke_idx = []
            for phase in range(number_of_phases):
                index_in_beat = np.floor(C[phase] * number_of_spokes_in_beat + np.arange(window_width))
                index_in_beat_cyclic = cyclic_spoke_axis[number_of_spokes_in_beat + index_in_beat.astype(int)]
                spoke_idx.append(spokes_in_beat[index_in_beat_cyclic])

            cardiac_phase[beat*window_width:(beat+1)*window_width, venc, :] = np.array(spoke_idx).T

    # data_ should be [1, col, lin, cha, 1, venc, phs]
    data_ = data[:, :, cardiac_phase, :]
    data_ = np.transpose(atleast_nd(data_, 7), (0, 2, 3, 6, 1, 4, 5))

    # traj_ should be [3, col, lin, 1, 1, venc, phs]
    traj_ = traj[:, :, cardiac_phase]
    traj_ = np.transpose(atleast_nd(traj_, 7), (2, 3, 4, 0, 1, 5, 6))

    # Reconstruct using BART (CG-SENSE with L2-regularization in image space)
    img = bart(1, f'pics -d5 -R Q:0.02', data_, sens, t=traj_)
    img = resize(img,(ncol//2, ncol//2, 1, 1, 1, number_of_vencs, number_of_phases))
    img_out = np.zeros((ncol//2, ncol//2, number_of_phases*2), dtype=complex)

    # Shared velocity encoding
    for ii in range(img.shape[-1]):
        z1 = img[..., 0, ii]
        z2 = img[..., 1, ii]
        if ii+1 < img.shape[-1]:
            z3 = img[..., 0, ii+1]
        else:
            z3 = img[..., 0, 0]

        # From Bernstein MA et. al. Reconstructions of phase contrast, phased array multicoil data.
        c1 = np.arctan2(np.imag(z1*np.conj(z2)), np.real(z1*np.conj(z2)))
        c2 = np.arctan2(np.imag(z3*np.conj(z2)), np.real(z3*np.conj(z2)))

        m1 = np.abs(z1) + np.abs(z2)
        m2 = np.abs(z3) + np.abs(z2)

        img_out[..., ii*2] = np.squeeze(m1 * np.exp( 1j * c1 ))
        img_out[..., ii*2+1] = np.squeeze(m2 * np.exp( 1j * c2 ))

    # From Gadgetron: mri_core/ScaleGadget.cpp and mri_core/FloatToFixedPointGadget.cpp
    # Magnitude images: Values above 4095 will be clamped.
    max_intensity_value = 4095
    min_intensity_value = 0
    intensity_offset_value = 2048

    data_mag = np.abs(img_out)
    data_mag *= max_intensity_value / data_mag.max()
    data_mag *= data_mag.max() / np.percentile(data_mag, 98.5)
    data_mag = np.around(data_mag)
    data_mag = np.clip(data_mag, min_intensity_value, max_intensity_value)
    data_mag = data_mag.astype(np.int16)

    rawHead = group[0].getHead()

    imagesOut = []

    for phs in range(data_mag.shape[-1]):
        tmpImg = ismrmrd.Image.from_array(data_mag[..., phs], transpose=False)

        # Set the header information
        tmpImg.setHead(update_img_header_from_raw(tmpImg.getHead(), rawHead))
        tmpImg.field_of_view = (ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.x),
                                ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.y),
                                ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.z))

        tmpImg.phase = phs
        tmpImg.image_index = phs + image_index_offset + 1 # Siemens begins image numbering at 1
        tmpImg.image_series_index = 1
        tmpImg.image_type = ismrmrd.IMTYPE_MAGNITUDE
        tmpImg.physiology_time_stamp[0] = ctypes.c_uint(triggerTime[phs])

        # Set ISMRMRD Meta Attributes
        tmpMeta = ismrmrd.Meta()
        tmpMeta['DataRole']               = 'Image'
        tmpMeta['ImageProcessingHistory'] = ['FIRE', 'SWIGPC']
        tmpMeta['WindowCenter']           = f'{intensity_offset_value}'
        tmpMeta['WindowWidth']            = f'{max_intensity_value}'
        tmpMeta['ImageComments']          = f'RR {RR_avg} +/- {RR_std}; {number_of_beats} heartbeats'

        # Add image orientation directions to MetaAttributes if not already present
        if tmpMeta.get('ImageRowDir') is None:
            tmpMeta['ImageRowDir'] = ["{:.18f}".format(rawHead.read_dir[0]), "{:.18f}".format(rawHead.read_dir[1]), "{:.18f}".format(rawHead.read_dir[2])]

        if tmpMeta.get('ImageColumnDir') is None:
            tmpMeta['ImageColumnDir'] = ["{:.18f}".format(rawHead.phase_dir[0]), "{:.18f}".format(rawHead.phase_dir[1]), "{:.18f}".format(rawHead.phase_dir[2])]

        xml = tmpMeta.serialize()
        tmpImg.attribute_string = xml
        imagesOut.append(tmpImg)

    # From Gadgetron: mri_core/FloatToFixedPointGadget.cpp
    # Phase: -pi will be 0, +pi will be 4095.
    data_phase = np.angle(img_out)
    data_phase *= intensity_offset_value / np.pi
    data_phase += intensity_offset_value
    data_phase = np.around(data_phase)
    data_phase = np.clip(data_phase, min_intensity_value, max_intensity_value)
    data_phase = data_phase.astype(np.int16)

    for phs in range(data_phase.shape[-1]):
        tmpImg = ismrmrd.Image.from_array(data_phase[..., phs], transpose=False)

        # Set the header information
        tmpImg.setHead(update_img_header_from_raw(tmpImg.getHead(), rawHead))
        tmpImg.field_of_view = (ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.x),
                                ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.y),
                                ctypes.c_float(metadata.encoding[0].reconSpace.fieldOfView_mm.z))

        tmpImg.phase = phs
        tmpImg.image_index = phs + 1
        tmpImg.image_series_index = 2
        tmpImg.image_type = ismrmrd.IMTYPE_PHASE
        tmpImg.physiology_time_stamp[0] = ctypes.c_uint(triggerTime[phs])

        # Set ISMRMRD Meta Attributes
        tmpMeta = ismrmrd.Meta()
        tmpMeta['DataRole']               = 'Phase'
        tmpMeta['ImageProcessingHistory'] = ['FIRE', 'SWIGPC']
        tmpMeta['WindowCenter']           = f'{intensity_offset_value}'
        tmpMeta['WindowWidth']            = f'{max_intensity_value}'
        tmpMeta['ImageComments']          = f'RR {RR_avg} +/- {RR_std}; {number_of_beats} heartbeats'
        tmpMeta['Keep_image_orientation'] = 1

        # Add image orientation directions to MetaAttributes if not already present
        if tmpMeta.get('ImageRowDir') is None:
            tmpMeta['ImageRowDir'] = ["{:.18f}".format(rawHead.read_dir[0]), "{:.18f}".format(rawHead.read_dir[1]), "{:.18f}".format(rawHead.read_dir[2])]

        if tmpMeta.get('ImageColumnDir') is None:
            tmpMeta['ImageColumnDir'] = ["{:.18f}".format(rawHead.phase_dir[0]), "{:.18f}".format(rawHead.phase_dir[1]), "{:.18f}".format(rawHead.phase_dir[2])]

        xml = tmpMeta.serialize()
        tmpImg.attribute_string = xml
        imagesOut.append(tmpImg)

    return imagesOut
