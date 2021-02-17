#!/bin/bash

old_version="1.1.5"
new_version="1.1.6"
mv relate_v${old_version}_MacOSX/ ./relate_v${new_version}_MacOSX/
mv relate_v${old_version}_x86_64_dynamic/ ./relate_v${new_version}_x86_64_dynamic/
mv relate_v${old_version}_x86_64_static/ ./relate_v${new_version}_x86_64_static/
