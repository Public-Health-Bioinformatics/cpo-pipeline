import os
import datetime
import drmaa
from structlog import get_logger

from cpo_pipeline.logging import now

logger = get_logger()

def prepare_job(job, session):
    """
    """
    job_template = session.createJobTemplate()
    job_template.jobName = job['job_name']
    job_template.nativeSpecification = job['native_specification']
    job_template.jobEnvironment = os.environ
    job_template.workingDirectory = os.getcwd()
    job_template.remoteCommand = job['remote_command']
    job_template.args = job['args']
    # job_template.joinFiles = True
    try:
        job_template.outputPath = ':' + job['output_path']
    except KeyError:
        pass
    try:
        job_template.errorPath = ':' + job['error_path']
    except KeyError:
        pass
    
    return job_template

def run_jobs(jobs):
    """
    """
    with drmaa.Session() as session:
        running_jobs = []
        for job in jobs:
            prepared_job = prepare_job(job, session)
            job_id = session.runJob(prepared_job)
            job_name = prepared_job.jobName
            logger.info(
                "job_submitted",
                timestamp=str(now()),
                job_name=job_name,
                job_id=job_id,
            )
            running_jobs.append({"id": job_id, "name": job_name})
        session.synchronize([x['id'] for x in running_jobs], drmaa.Session.TIMEOUT_WAIT_FOREVER, False)
        
        for job in running_jobs:
            job_info = session.wait(job["id"], drmaa.Session.TIMEOUT_WAIT_FOREVER)
            resource_usage = job_info.resourceUsage
            float_fields = [
                "io",
                "iow",
                "mem",
                "cpu",
                "vmem",
                "maxvmem",
                "priority",
                "ru_wallclock",
                "ru_utime",
                "ru_stime",
                "ru_maxrss",
                "ru_ixrss",
                "ru_ismrss",
                "ru_idrss",
                "ru_isrss",
                "ru_minflt",
                "ru_majflt",
                "ru_nswap",
                "ru_inblock",
                "ru_oublock",
                "ru_msgsnd",
                "ru_msgrcv",
                "ru_nsignals",
                "ru_nvcsw",
                "ru_nivcsw",
                "acct_cpu",
                "acct_mem",
                "acct_io",
                "acct_iow",
                "acct_maxvmem",
            ]
            for float_field in float_fields:
                resource_usage[float_field] = float(resource_usage[float_field])
            int_fields = [
                "exit_status"
            ]
            for int_field in int_fields:
                resource_usage[int_field] = int(float(resource_usage[int_field]))
            # Convert unix epoch timestamps to ISO8601 (YYYY-MM-DDTHH:mm:ss+tz)
            time_fields = [
                "submission_time",
                "start_time",
                "end_time",
            ]
            for time_field in time_fields:
                unix_timestamp = resource_usage[time_field]
                iso8601_timestamp = str(datetime.datetime.fromtimestamp(
                    int(float(unix_timestamp)), datetime.timezone.utc
                ).isoformat())
                resource_usage[time_field] = iso8601_timestamp
            logger.info(
                "job_completed",
                timestamp=str(now()),
                job_id=job["id"],
                job_name=job["name"],
                resource_usage=job_info.resourceUsage,
                exit_status=job_info.exitStatus,
            )
