'''Function to provide update bar on terminal
during loop processes. Somewhat shamelessly
snatched form the internet but luckily Python
is open source.

Author: Not Simon Zieleniewski

Last updated: 25-11-15

'''

import sys


def update_progress(progress):
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
##    if not isinstance(progress, float):
##        progress = 0
##        status = "error: progress var must be float\r\n"
##    if progress < 0:
##        progress = 0
##        status = "Halt...\r\n"
##    if progress >= 1:
##        progress = 1
##        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "="*block + " "*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
