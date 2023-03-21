import argparse
import json
import yaml
from augur.utils import write_json

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--auspice-config-path', type=str)
    parser.add_argument('--snake-config-path', type=str)
    parser.add_argument('--output-config-path', type=str)

    args = parser.parse_args()

    auspice_config = json.load(open(args.auspice_config_path, "r"))
    snake_config = yaml.safe_load(open(args.snake_config_path, "r"))

    if 'polyclonal_serum_models' in snake_config.keys():
        for exp_label, exp_info in snake_config['polyclonal_serum_models'].items():

            if 'concentrations' in exp_info:
                for c in exp_info['concentrations'].split(","):
                    auspice_config['colorings'].append({
                        "key" : f"{exp_label}_prob_escape_c_{c}",
                        "title" : f"{exp_label}_prob_escape_c_{c}",
                        "type": "continuous",
                        "scale": [[0.0, "#2B83BA"], [0.5, "#FBC93D"], [1.0, "#D7191C"]],
                        "legend": [
                            {"value": 0.001, "display": "0.00", "bounds": [0.000, 0.125]},
                            {"value": 0.25, "display": "0.25", "bounds": [0.125, 0.375]},
                            {"value": 0.50, "display": "0.50", "bounds": [0.375, 0.625]},
                            {"value": 0.75, "display": "0.75", "bounds": [0.625, 0.875]},
                            {"value": 1.00, "display": "1.00", "bounds": [0.875, 1.000]}
                        ]
                    })

                    
            if 'icXX' in exp_info:
                col = f"ic{int(exp_info['icXX']*100)}"
                auspice_config['colorings'].append({
                    "key" : f"{exp_label}_{col}",
                    "title" : f"{exp_label}_{col}",
                    "type": "continuous",
                    # TODO are these scales still going to be correct?
                    #"scale": [[0.0, "#2B83BA"], [0.5, "#FBC93D"], [1.0, "#D7191C"]],
                    #"legend": [
                    #    {"value": 0.001, "display": "0.00", "bounds": [0.000, 0.125]},
                    #    {"value": 0.25, "display": "0.25", "bounds": [0.125, 0.375]},
                    #    {"value": 0.50, "display": "0.50", "bounds": [0.375, 0.625]},
                    #    {"value": 0.75, "display": "0.75", "bounds": [0.625, 0.875]},
                    #    {"value": 1.00, "display": "1.00", "bounds": [0.875, 1.000]}
                    #]
                })


    if 'escape_models' in snake_config.keys():
        for exp_label, exp_info in snake_config['escape_models'].items():
            auspice_config['colorings'].append({
                "key" : f"{exp_label}_additive_escape",
                "title" : f"{exp_label}_additive_escape",
                "type": "continuous",
                "scale": [[0.0, "#2B83BA"], [0.5, "#FBC93D"], [1.0, "#D7191C"]],
                "legend": [
                    {"value": 0.001, "display": "0.00", "bounds": [0.000, 0.125]},
                    {"value": 0.25, "display": "0.25", "bounds": [0.125, 0.375]},
                    {"value": 0.50, "display": "0.50", "bounds": [0.375, 0.625]},
                    {"value": 0.75, "display": "0.75", "bounds": [0.625, 0.875]},
                    {"value": 1.00, "display": "1.00", "bounds": [0.875, 1.000]}
                ]
            })

    write_json(auspice_config, args.output_config_path, include_version=False)
