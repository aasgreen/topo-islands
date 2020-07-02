#!/bin/bash

docker run --memory-swap -1 -it -v "${PWD}:/work/" aasgreen/topo-simulations:1.4 ./call_soliton.sh
