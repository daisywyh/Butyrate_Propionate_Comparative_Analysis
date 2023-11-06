import sys
o_file = sys.argv[1]
n_file = sys.argv[2]

with open(o_file, 'r') as f_read:
    with open(n_file, 'w') as f_write:
        i = 1
        for line in f_read:
            if (i % 2) != 0:
                n_line = ">" + line
                f_write.write(n_line)
            elif (i % 2) == 0:
                f_write.write(line)
            i += 1
