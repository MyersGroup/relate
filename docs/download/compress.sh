#!/bin/bash

version="1.2.1"

tar -czvf ./relate_v${version}_MacOSX_Intel.tgz relate_v${version}_MacOSX_Intel
tar -czvf ./relate_v${version}_MacOSX_M.tgz relate_v${version}_MacOSX_M
tar -czvf ./relate_v${version}_x86_64_dynamic.tgz relate_v${version}_x86_64_dynamic
tar -czvf ./relate_v${version}_x86_64_static.tgz relate_v${version}_x86_64_static

