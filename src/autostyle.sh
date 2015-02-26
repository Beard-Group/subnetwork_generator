#!/bin/bash
# Once astyle is installed (check aliases), run "source auto_style.sh" on command line

astyle --style=stroustrup --convert-tabs --max-code-length=80 \
    main.cu
