#!/bin/bash

cp /home/tomo/Projet_tomo/Tomo_config/recon.txt /ramdisk/ACQUIS/ && tomo_pretraitement -i /ramdisk/ACQUIS/ -o /home/tomo/Projet_tomo/ && tomo_reconstruction
#reconstruction




