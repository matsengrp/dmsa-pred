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
                    auspice_config['colorings'].insert(0, {
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

                    
            if 'icxx' in exp_info:
                col = f"ic{int(exp_info['icxx']*100)}_log_fold_change"
                auspice_config['colorings'].insert(0, {
                    "key" : f"{exp_label}_{col}",
                    "title" : f"{exp_label}_{col}",
                    "type": "continuous",
                    "scale": [
                        [-2.0, "#E1C1F1"], [0.0, "#2B83BA"], [2, "#FBC93D"], [4, "#D7191C"]
                    ],
                    "legend": [
                        {"value": -1.0, "display": "-1.0", "bounds": [-100.0, -0.5]},
                        {"value": 0.0, "display": "0.0", "bounds": [-0.5, 0.5]},
                        {"value": 1.0, "display": "1.0", "bounds": [0.5, 1.5]},
                        {"value": 2.0, "display": "2.0", "bounds": [1.5, 2.5]},
                        {"value": 3.0, "display": "3.0", "bounds": [2.5, 3.5]},
                        {"value": 4.0, "display": "4.0", "bounds": [3.5, 100.0]},
                    ]
                })


    if 'escape_models' in snake_config.keys():
        for exp_label, exp_info in snake_config['escape_models'].items():
            auspice_config['colorings'].insert(0, {
                "key" : f"{exp_label}",
                "title" : f"{exp_label}",
                "type": "continuous",
                    "scale": [
                        [0.0, "#2B83BA"], [0.5, "#FBC93D"], [1.0, "#D7191C"]
                    ],
                    "legend": [
                        {"value": 0.0, "display": "0.0", "bounds": [-100, 0.125]},
                        {"value": 0.25, "display": "0.25", "bounds": [0.125, 0.375]},
                        {"value": 0.5, "display": "0.5", "bounds": [0.375, 0.625]},
                        {"value": 0.75, "display": "0.75", "bounds": [0.625, 0.875]},
                        {"value": 1.0, "display": "1.0", "bounds": [0.875, 100]}
                    ]
            })

    write_json(auspice_config, args.output_config_path, include_version=False)
