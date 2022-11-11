# SWIG

FIRE implementation of the reconstruction described in Fyrdahl et. al. Magn Reson Med. 2020 Apr;83(4):1310-1321.

## Prerequisites
---
:whale: To run this code, you need to install [Docker](https://docs.docker.com/get-docker/).


## How-to
---
:zap: To quickly get started, you can run a reconstruction on the example dataset with the command
```
docker compose up
```
This will build the server (from local code), pull an example dataset from Amazon S3, and run the reconstruction. The resulting images can be found in the file out.h5 in the root folder.


:hand: If you prefer a more hands-on approach, you can pull the server from Dockerhub using
```
docker pull fyrdahl/swig
````
or build the server from code using
```
docker build --tag=swig_dev -f docker/Dockerfile .
```

To run the server in detached mode, with the current working directory $(pwd) mounted to /tmp, run
```
docker run -p=9002:9002 --rm --detach --volume=$(pwd):/tmp --name swig_server swig_dev
```

To get access to the the client script, you can execute an interactive shell in the already running container
```
docker exec -ti swig_server /bin/bash
```

To reconstruct the radial through-plane phase-contrast dataset, run the following command in the interactive client session
```
python3 client.py -p 9002 -o /tmp/mitral_tp_recon.h5 -c swig /data/mitral_tp.h5
```

Likewise, to reconstruct the radial in-plane phase-contrast dataset, run the following command
```
python3 client.py -p 9002 -o /tmp/mitral_ip_recon.h5 -c swig /data/mitral_ip.h5
```


## References
---
:bulb: If you found this code useful, consider citing:
>``` Fyrdahl A, Ramos JG, Eriksson MJ, Caidahl K, Ugander M, Sigfridsson A. Sector-wise golden-angle phase contrast with high temporal resolution for evaluation of left ventricular diastolic dysfunction. Magn Reson Med. 2020 Apr;83(4):1310-1321. doi: 10.1002/mrm.28018```

Files for import into common reference managers are available [here](https://github.com/fyrdahl/SWIG/tree/main/ref)

List of other references relevant to this code;
* Rosenzweig S, Holme HCM, Uecker M. Simple auto-calibrated gradient delay estimation from few spokes using Radial Intersections (RING). Magn Reson Med. 2019; 81:1898–1906.
* Inati SJ, Hansen MS, Kellman P. A Solution to the Phase Problem in Adaptive Coil Combination. Annual Meeting ISMRM, Salt Lake City 2013, In Proc. Intl. Soc. Mag. Reson. Med. 21:2672
* Bernstein MA, Grgic M, Brosnan TJ, Pelc NJ, Reconstructions of phase contrast, phased array multicoil data. Magn. Reson. Med. 1994; 32:330-334.
* Uecker M, Ong F, Tamir JI, Bahri D, Virtue P, Cheng JY, Zhang T, Lustig M. Berkeley Advanced Reconstruction Toolbox, Annual Meeting ISMRM, Toronto 2015, In Proc. Intl. Soc. Mag. Reson. Med. 23:2486
