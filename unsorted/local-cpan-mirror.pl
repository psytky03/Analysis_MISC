Forked from: https://gist.github.com/reyjrar/1372269

#=======================================
# Part 1 is Setting up the Mirror Server

# Install CPAN::Mini
$ curl -L http://cpanmin.us | perl - --sudo CPAN::Mini

# Select a CPAN Mirror URL from http://mirrors.cpan.org/
#   - We'll use http://cpan.pair.com
# Pick a directory to mirror to, I'll use /var/www/cpan

# Run minicpan, for more details see: https://metacpan.org/module/minicpan
$ minicpan -l /var/www/cpan -r http://cpan.pair.com


#=======================================================
# Part 2 is creating a CPAN/Config.pm to use your mirror
#
#  We'll assume we mirrored off of dev.example.com with 
#  a simple subdirectory called cpan.  Our URL is:
#    http://dev.example.com/cpan

# Now, on a CLIENT machine:
$ sudo cpan
cpan> o conf urllist "http://dev.example.com/cpan"
cpan> o conf prerequisites_policy "follow"
cpan> o conf commit
cpan> exit

# Now that client can install modules from your local mirror
#  - Example, installing Mouse:
$ cpan Mouse

#==========================
# Part 3 Super Bonus Round
#
# Oh, you are a clever bastard aren't you.  You want to carry the CPAN in your pocket!
#  Well played, good sir.
#
# Insert a Flash Drive with 4GB or Greater Capacity!

$ minicpan -l /Volumes/SuperDuperFlashyDrive/cpan -r http://cpan.pair.com

# Carry your fancy Thumb Drive to your offline computer.
# Plug it in.
# It mounts to /media/SuperDuperFlashyDrive
# Get CRAZY!
$ cpan
cpan> o conf urllist "file:///media/SuperDuperFlashyDrive/cpan"
cpan> o conf commit
cpan> install Mouse

# There you have it.. cpanminus, minicpan, and a portable CPAN!
