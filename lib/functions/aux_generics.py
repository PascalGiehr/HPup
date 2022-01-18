# coding=utf-8

"""
Module containing generic functions that can be part of any
pipeline
"""

import os as os
import re as re
import subprocess as sp
import fnmatch as fnm
import time as time

from ruffus import JobSignalledBreak


class SafeReplace(dict):
    def __missing__(self, key):
        return '{' + key + '}'


def _normalize_job_output(output):
    """
    :param output:
    :return:
    """
    if isinstance(output, list):
        if len(output) == 0:
            normout = ''
        elif isinstance(output[0], bytes):
            output = list(map(lambda x: x.decode('utf-8'), output))
            normout = '\n'.join(output).strip()
        else:
            assert isinstance(output[0], str), 'Unexpected type of Job output (in list): {}'.format(type(output[0]))
            normout = '\n'.join(output).strip()
    elif isinstance(output, bytes):
        normout = output.decode('utf-8').strip()
    else:
        assert isinstance(output, str), 'Unexpected type of Job output: {}'.format(type(output))
        normout = output.strip()
    return normout


def check_job(out, err):
    """
    :param out:
    :param err:
    :return:
    """
    error_notes = ['error', 'fail', 'failed', 'failure', 'segfault', 'abort']
    out = _normalize_job_output(out)
    err = _normalize_job_output(err)
    if err:
        if any([msg in err.lower() for msg in error_notes]):
            try:
                from ruffus import JobSignalledBreak
                raise JobSignalledBreak(err)
            except ImportError:
                raise RuntimeError(err)
    return out, err


def assert_exists(outputfile):
    """
    :param outputfile:
    :return:
    """
    assert os.path.isfile(outputfile), 'Given outputfile {} is not a regular file or link'.format(outputfile)
    return outputfile


def syscall_raw(cmd, syscall):
    """
    :param cmd:
    :param syscall:
    :return:
    """
    out, err = syscall(cmd)
    out, err = check_job(out, err)
    return


def syscall_in_out(inputfile, outputfile, cmd, syscall, posrep=False):
    """
    :param inputfile:
    :param outputfile:
    :param cmd:
    :param syscall:
    :param posrep: use positional replacement for input and output when formatting command
    :return:
    """
#    print("Inputfile: "+str(inputfile))
    assert os.path.isfile(inputfile), 'Input path is not a file: {}'.format(inputfile)
    assert outputfile, 'Received no output file'
    if posrep:
        cmd = cmd.format(*(inputfile, outputfile))
    else:
        files = {'inputfile': inputfile, 'outputfile': outputfile}
        cmd = cmd.format(**files)
#    print("Command: "+cmd)
    out, err = syscall(cmd)
    out, err = check_job(out, err)
    assert os.path.isfile(outputfile), 'Output path is not a file: {} - job failed?'.format(outputfile)
    return outputfile

def syscall_ins_outs(inputfiles, outputfiles, cmd, syscall, posrep=False):
    """
    :param inputfiles:
    :param outputfiles:
    :param cmd:
    :param syscall:
    :param posrep: use positional replacement for input and output when formatting command
    :return:
    """
 
    if type(inputfiles) == type(""):
        inputfiles= [inputfiles]

    #print("!inputs!")
    #print(type(inputfiles))
    #print(inputfiles)
    #print(outputfiles)
    #print(cmd)
    
    assert all([os.path.isfile(f) for f in inputfiles]), 'Not all input paths are files: {}'.format(inputfiles)
    assert outputfiles, 'Received no output file'
    #print(cmd)
    if posrep:
        cmd = cmd.format(*(inputfiles, outputfiles))
    else:
        files = {'inputfile': inputfiles, 'outputfile': outputfiles}
        cmd = cmd.format(**files)
    #print("Cmd: "+cmd)
    out, err = syscall(cmd)
    #print("What is out?")
    #print(out)
    #print("What is err?")
    #print(err)
    out, err = check_job(out, err)
    assert all([os.path.isfile(f) for f in outputfiles]), 'Not all output paths are files: {} - job failed?'.format(outputfiles)

    return outputfiles

def syscall_in_pat(inputfile, outputpattern, outdir, cmd, syscall, posrep=False):
    """ System call for cases where a single input file is split
    into multiple output files (number determined at runtime), hence
    outputfiles represents a matching pattern rather than a filename

    :return: list of files
    """
    assert os.path.isfile(inputfile), 'Input path is not a file: {}'.format(inputfile)
    if posrep:
        cmd = cmd.format(*(inputfile,))
    else:
        cmd = cmd.format(**{'inputfile': inputfile})
    #print(cmd)
#     out, err = syscall(cmd)
#     out, err = check_job(out, err)
    outfiles = os.listdir(outdir)
    outfiles = [os.path.join(outdir, f) for f in outfiles]
    if isinstance(outputpattern, list) and len(outputpattern) == 1:
        outputpattern = str(outputpattern[0])
        outfiles = fnm.filter(outfiles, outputpattern)
    elif isinstance(outputpattern, str) and len(outputpattern) > 0:
        outfiles = fnm.filter(outfiles, outputpattern)
    else:
        pass
    assert len(outfiles) > 0, 'No output files produced by command {} with filter pattern {}'.format(cmd, outputpattern)
    return outfiles


def syscall_in_out_ref(inputfile, outputfile, cmd, reffile):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    assert os.path.isfile(inputfile), 'Input path is not a file: {}'.format(inputfile)
    assert outputfile, 'Received no output file'
    assert os.path.isfile(reffile), 'Reference path is not a file: {}'.format(outputfile)
    files = {'inputfile': inputfile, 'outputfile': outputfile, 'reffile': reffile}
    cmd = cmd.format(**files)
    out, err = syscall(cmd)
    out, err = check_job(out, err)
    assert os.path.isfile(outputfile), 'Output path is not a file: {} - job failed?'.format(outputfile)
    return outputfile


def syscall_ins_out(inputfiles, outputfile, cmd, syscall, posrep=False):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    if type(inputfiles) == type(""): 
        inputfiles= [inputfiles]
    #print("The input:")
    #print(inputfiles)    
    #print("command:")
    #print(cmd)
    #print("output:")
    #print(outputfile)   
    assert all([os.path.isfile(f) for f in inputfiles]), 'Not all input paths are files: {}'.format(inputfiles)
    assert outputfile, 'Received no output file'
    if posrep:
        cmd = cmd.format(*(outputfile,))
    else:
        cmd = cmd.format(**{'outputfile': outputfile})
    #print(cmd)
    out, err = syscall(cmd)
    out, err = check_job(out, err)
    assert os.path.isfile(outputfile), 'Output path is not a file: {} - job failed?'.format(outputfile)
    return outputfile

def syscall_ins_out_KN(inputfiles, outputfile, cmd, syscall, posrep=False, separator=" "):
    """
    :param inputfiles:
    :param outputfile:
    :return:
    """
    #if type(inputfiles) == type(""):
     #   inputfiles= [inputfiles]
    print("inputs")
    print(inputfiles)   
    assert all([os.path.isfile(f) for f in inputfiles]), 'Not all input paths are files: {}'.format(inputfiles)
    assert outputfile, 'Received no output file'
    if posrep:
        # change here input
        cmd = cmd.format(*(separator.join(inputfiles), outputfile))
    else:
        files = {'inputfile': separator.join(inputfiles), 'outputfile': outputfile}
        cmd = cmd.format(**files)
    out, err = syscall(cmd)
    out, err = check_job(out, err)
    #print("cmd")
    #print(cmd)
    assert os.path.isfile(outputfile), 'Output path is not a file: {} - job failed?'.format(outputfile)
    return outputfile

def syscall_ins_pat(inputfiles, outputpattern, outdir, cmd, syscall, posrep=False):
    """ System call for cases where a set of input files is split
    into multiple output files (number determined at runtime), hence
    outputfiles represents a matching pattern rather than a filename

    :return: list of files
    """
    
    assert all([os.path.isfile(f) for f in inputfiles]), 'Not all input paths are files: {}'.format(inputfiles)
    assert os.path.isdir(outdir), 'Output dir is not a folder: {}'.format(outdir)
    if posrep:
        cmd = cmd.format(*(str(inputfiles),))
    else:
        cmd = cmd.format(**{'inputfiles': str(inputfiles)})
    out, err = syscall(cmd)
    out, err = check_job(out, err)
    outfiles = os.listdir(outdir)
    outfiles = [os.path.join(outdir, f) for f in outfiles]
    if isinstance(outputpattern, list) and len(outputpattern) == 1:
        outputpattern = str(outputpattern[0])
        outfiles = fnm.filter(outfiles, outputpattern)
    elif isinstance(outputpattern, str) and len(outputpattern) > 0:
        outfiles = fnm.filter(outfiles, outputpattern)
    else:
        pass
    assert len(outfiles) > 0, 'No output files produced by command {} with filter pattern {}'.format(cmd, outputpattern)
    return outfiles

def syscall_ins_pat_KN(inputfiles, outputpattern, outdir, cmd, syscall, posrep=False):
    """ System call for cases where a set of input files is split
    into multiple output files (number determined at runtime), hence
    outputfiles represents a matching pattern rather than a filename

    :return: list of files
    """
    
    assert all([os.path.isfile(f) for f in inputfiles]), 'Not all input paths are files: {}'.format(inputfiles)
    assert os.path.isdir(outdir), 'Output dir is not a folder: {}'.format(outdir)
    if posrep:
        cmd = cmd.format(*inputfiles)
    else:
        cmd = cmd.format(**{'inputfiles': str(inputfiles)})
    out, err = syscall(cmd)
    out, err = check_job(out, err)
    outfiles = os.listdir(outdir)
    outfiles = [os.path.join(outdir, f) for f in outfiles]
    if isinstance(outputpattern, list) and len(outputpattern) == 1:
        outputpattern = str(outputpattern[0])
        outfiles = fnm.filter(outfiles, outputpattern)
    elif isinstance(outputpattern, str) and len(outputpattern) > 0:
        outfiles = fnm.filter(outfiles, outputpattern)
    else:
        pass
    outfiles.sort()
#     print(outfiles)
    
    assert len(outfiles) > 0, 'No output files produced by command {} with filter pattern {}'.format(cmd, outputpattern)
    return outfiles

def syscall(cmdline, log, lock, env=None):
    """
    :param cmdline:
    :param log:
    :param lock:
    :param env:
    :return:
    """
    try:
        _ = sp.check_call(cmdline, shell=True, env=env)
    except sp.CalledProcessError as spe:
        msg = 'System call {} returned error code: {}'.format(spe.cmd, spe.returncode)
        with lock:
            log.error(msg)
        raise JobSignalledBreak('Error: {} failed'.format(cmdline))
    else:
        return


def collect_files_from_tree(root, filter=None):
    """ Traverse tree and collect full paths to all files. This does not
    follow links

    :param root: start traversal here
     :type: str
    :param filter: filter list for files matching filter, e.g. *.txt
     :type: str
    :return: a list containing full paths to all files in the subtree
     :rtype: list of str
    """
    assert os.path.isdir(root), 'Given root is not a folder: {}'.format(root)
    filelist = []
    for dp, dn, fn in os.walk(root):
        if fn:
            for f in fn:
                filelist.append(os.path.join(dp, f))

    if filter:
        filelist = fnm.filter(filelist, filter)
    return filelist


def link_input_data(topdir, outdir, patterns):
    """ Symlink all files matching any of patterns in outdir. Do
    not overwrite/recreate already existing links
    :param topdir:
    :param outdir:
    :param patterns:
    :return:
    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    cpats = [re.compile(p) for p in patterns]
    links = []
    for root, dirs, files in os.walk(topdir):
        if files:
            for f in files:
                if any([p.search(f) is not None for p in cpats]):
                    trg = os.path.join(outdir, f)
                    if os.path.isfile(trg) or os.path.islink(trg):
                        links.append(trg)
                        continue
                    src = os.path.join(root, f)
                    os.symlink(src, trg)
                    links.append(trg)
    assert len(links) > 0, 'No files linked for patterns {} in tree {}'.format(patterns, topdir)
    return links


def check_tool_versions(version_infos, config, sci_obj, outputfile):
    """
    :param inputfiles:
    :param outputfile:
    :param config:
    :param sci_obj:
    :return:
    """
    if sci_obj is None:
        with open(outputfile, 'w') as outfile:
            pass
        return
    pathconf = {'workdir': os.path.dirname(outputfile)}
    infos = []
    it = iter(version_infos)
    for entry in it:
        env = entry
        call = next(it)
        check = next(it)
        toolname = call[0].split('_')[0]
        assert env[0].startswith(toolname) and check[0].startswith(toolname),\
            'Tool name mismatch during version check'
        exenv = dict(config.items(env[1]))
        sci_obj.set_config_env(pathconf, exenv)
        syscall = sci_obj.local_job()
        out, err = syscall(call[1])
        if 'ERROR' in err or 'error' in out.lower():
            pperr = err.replace('\t', ' ')
            pperr = pperr.replace('\n', ' - ')
            infos.append('{}\t{}'.format(toolname, 'non-zero exit: {}'.format(pperr)))
        elif check[1] in out or check[1] in err:
            infos.append('{}\t{}'.format(toolname, check[1]))
        else:
            ppout = out.replace('\t', ' ')
            ppout = ppout.replace('\n', ' - ')
            pperr = err.replace('\t', ' ')
            pperr = pperr.replace('\n', ' - ')
            total = (ppout + ' ' + pperr).strip()
            infos.append('{}\t{}: {}'.format(toolname, 'no_match', total))
    with open(outputfile, 'w') as outfile:
        _ = outfile.write('\n'.join(infos))
        _ = outfile.write('\n')
    return outputfile


def update_time_sheet(timesheet):
    """
    :param timesheet:
    :return:
    """
    with open(timesheet, 'a') as outfile:
        _ = outfile.write(str(int(time.time())) + '\n')
        _ = outfile.write(time.ctime() + '\n')
    return
