#!/usr/bin/awk -f
BEGIN { start = 0; useThreshold = 0; threshold = 0.0 };
/------ inclusion threshold ------/ { if (useThreshold) { next } else { exit } };
/^$/ && start > 0 { exit }; # exit on empty line
/Scores for complete hits/ { start = 1; next }; # found the start of the table, need to skip a couple of lines
start > 0 && start < 3 { start++; next }; # skip over lines till we reach the third
start == 3 { # amongst the data
        if (useThreshold) { # using a custom threshold
                if ($1 >= threshold) next
        }
        if ($5 < $6) {
                printf("%s\t%d\t%d\t%s\tforward\t+\n", $4, $5, $6, $4)
        }
        else {
                printf("%s\t%d\t%d\t%s\treverse\t-\n", $4, $6, $5, $4)
        }
};
