#!/bin/bash

version="1.0.17"

tar -czvf ./relate_v${version}_MacOSX.tgz relate_v${version}_MacOSX
tar -czvf ./relate_v${version}_x86_64_dynamic.tgz relate_v${version}_x86_64_dynamic
tar -czvf ./relate_v${version}_x86_64_static.tgz relate_v${version}_x86_64_static

