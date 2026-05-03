#!/bin/sh
ldpath="${LD_LIBRARY_PATH:-$(</etc/ld.so.conf)}"
notfound=1
for libdir in ${ldpath//:/ }; do
        (test -f "$libdir/${1}" && echo "$_") && notfound=0
done
[ "$notfound" -eq 0 ]
