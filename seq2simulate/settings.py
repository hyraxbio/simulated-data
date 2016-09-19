import os
import json


def write_settings(folder, settings, samples):
    settings["sample_details"] = samples

    with open(os.path.join(folder, "settings.json"), 'w') as f:
        f.write(json.dumps(settings, indent=2)) 

