#!/usr/bin/env python

import os
import argparse
import configparser
import csv
import subprocess
import asyncio

import cpo_pipeline

def prepare_job(job, session):
    job_template = session.createJobTemplate()
    job_template.jobName = job['job_name']
    job_template.nativeSpecification = job['native_specification']
    job_template.jobEnvironment = os.environ
    job_template.workingDirectory = os.getcwd()
    job_template.remoteCommand = job['remote_command']
    job_template.args = job['args']
    job_template.joinFiles = True
    
    return job_template

async def run_command(*args):
    """Run command in subprocess
    
    Example from:
        http://asyncio.readthedocs.io/en/latest/subprocess.html
    """
    # Create subprocess
    process = await asyncio.create_subprocess_exec(
        *args,
        # stdout must be a pipe to be accessible as process.stdout
        stdout=asyncio.subprocess.PIPE
    )

    # Status
    print('Started:', args, '(pid = ' + str(process.pid) + ')')

    # Wait for the subprocess to finish
    stdout, stderr = await process.communicate()

    # Progress
    if process.returncode == 0:
        print('Done:', args, '(pid = ' + str(process.pid) + ')')
    else:
        print('Failed:', args, '(pid = ' + str(process.pid) + ')')

    # Result
    result = stdout.decode().strip()

    # Return stdout
    return result

def make_chunks(l, n):
    """Yield successive n-sized chunks from l.

    Note:
        Taken from https://stackoverflow.com/a/312464
    """
    # Assume Python 3
    for i in range(0, len(l), n):
        yield l[i:i + n]

def run_asyncio_commands(tasks, max_concurrent_tasks=0):
    """Run tasks asynchronously using asyncio and return results

    If max_concurrent_tasks are set to 0, no limit is applied.
    """

    all_results = []

    if max_concurrent_tasks == 0:
        chunks = [tasks]
    else:
        chunks = make_chunks(l=tasks, n=max_concurrent_tasks)

    for tasks_in_chunk in chunks:
        loop = asyncio.get_event_loop()

        commands = asyncio.gather(*tasks_in_chunk)  # Unpack list using *
        results = loop.run_until_complete(commands)
        all_results += results
        loop.close()
    return all_results


def main(args):
    """
    """
    
    config = configparser.ConfigParser()
    config.read(args.config_file)

    sample_ids = []
    with open(args.input_file) as input_file:
        fieldnames = ['sample_id', 'reads1_fastq', 'reads2_fastq']
        reader = csv.DictReader(
            (row for row in input_file if not row.startswith('#')),
            delimiter='\t',
            fieldnames=fieldnames
        )
        tasks = []
        for row in reader:
            sample_ids.append(row['sample_id'])
            command_line = [
                'cpo-pipeline',
                'assembly',
                '--ID', row['sample_id'],
                '--R1', row['reads1_fastq'],
                '--R2', row['reads2_fastq'],
                '--outdir', args.outdir,
            ]
            tasks.append(run_command(*command_line))
        results = run_asyncio_commands(tasks, max_concurrent_tasks=2)
        print('Results:', results)
        
    tasks = []
    for sample_id in sample_ids:
        command_line = [
            'cpo-pipeline',
            'typing',
            '--ID', sample_id,
            '--assembly', "/".join(args.outdir, sample_id, 'assembly', 'contigs.fa'),
            '--outdir', args.outdir,
        ]
        tasks.append(run_command(*command_line))
    results = run_asyncio_commands(tasks, max_concurrent_tasks=2)
    print('Results:', results)
        
if __name__ == '__main__':
    script_name = os.path.basename(os.path.realpath(sys.argv[0]))
    parser = argparse.ArgumentParser(prog=script_name, description='')
    parser.add_argument('-i', '--input',  dest='input_file',
                        help='Multi-Sample Input File', required=True)
    parser.add_argument('-o', '--outdir', dest='outdir',
                        help='Output Directory', required=True)
    parser.add_argument('-c', '--config', dest='config_file',
                        help='Config File', required=False)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {}'.format(cpo_pipeline.__version__))
    args = parser.parse_args()
    main(args)
