#!/bin/bash

bash codeml_loop.sh acceptance/ir_seq ir

bash codeml_loop.sh location/obp_seq obp
bash codeml_loop.sh location/or_seq or

bash codeml_loop.sh detox/abc_seq abc
bash codeml_loop.sh detox/cyp_seq cyp
bash codeml_loop.sh detox/est_seq est
bash codeml_loop.sh detox/gst_seq gst
bash codeml_loop.sh detox/udp_seq udp
