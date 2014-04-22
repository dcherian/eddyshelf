#!/usr/bin/python3


def make_montage(basedir, depths):
    """ makes a montage of passive tracer animation from runs.animate_pt
    run with different depths

    Arguments:
        basedir - basedir to which depths are appended i.e., runew-03-pt-z-
        depths  - depths at which stuff has been outputted

    Returns
         none

    Deepak Cherian - 23/01/2014

    """

    import subprocess
    import glob
    import os

    # first find number of files
    flist = glob.glob(basedir + str(depths[0]) + '/*.png')
    N = len(flist)

    print(depths)

    outdir = 'temp_pt'
    outfmt_av = './' + outdir + '/output_%06d.png'
    outfmt_mo = './' + outdir + '/output_{0:06d}.png'
    # get runname by splitting by '/'
    # and partitioning at the -z introduced by runs.animate_pt
    outname = basedir.split('/')[-1].rpartition('-z')[0] + '.mp4'

    # avconv options
    frameRate = 5
    bitRate = 25000000
    # avconvArgs = ''

    # make temp dir
    try:
        os.mkdir(outdir)
    except os.FileExistsError:
        subprocess.call(['rm', '-rf', outdir])
        os.mkdir(outdir)

    for ii in range(1, N+1):  # range(1,N):
        print('Processing image ' + str(ii) + '/' + str(N))
        fname = '/mm_frame_{0:06d}.png'.format(ii)

        # arguments for montage command
        argument = 'montage '
        for jj in depths:
            argument += basedir + str(jj) + fname + ' '
        argument += '-geometry 1600x900 ' + outfmt_mo.format(ii)

        # call the montage command for each set of images
        subprocess.call(argument.split())

    # all output images have been created
    # now execute avconv command
    avconv = ('avconv -r {0} -f image2 -i {1} -q:v 1 -g 1 -b:v {2} {3}'.
              format(frameRate, outfmt_av, bitRate, outname))
    print(avconv)
    subprocess.call(avconv.split())

if __name__ == '__main__':
    import sys

    make_montage(sys.argv[1], sys.argv[2:])
