#!/bin/bash

cd dataPublic
mkdir -p UXM
cd UXM

url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05580-6/MediaObjects/41586_2022_5580_MOESM5_ESM.zip"
curl "${url}"
gunzip "41586_2022_5580_MOESM5_ESM.zip"