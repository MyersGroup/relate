#!/bin/bash

old_version="1.1.9"
new_version="1.2.0"
mv relate_v${old_version}_MacOSX_Intel/ ./relate_v${new_version}_MacOSX_Intel/
mv relate_v${old_version}_MacOSX_M/ ./relate_v${new_version}_MacOSX_M/
mv relate_v${old_version}_x86_64_dynamic/ ./relate_v${new_version}_x86_64_dynamic/
mv relate_v${old_version}_x86_64_static/ ./relate_v${new_version}_x86_64_static/
