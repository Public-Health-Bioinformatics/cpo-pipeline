import datetime

def now():
    return datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).isoformat()
