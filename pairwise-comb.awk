#!/usr/bin/awk -f

#Cartesian product of records

{
    file = FILENAME
    while ((getline line <file) > 0)
        print $0, line
    close(file)
}
