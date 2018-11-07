cpp -P -fpreprocessed | tr -d '[:space:]' | md5sum | cut -d ' ' -f 1
