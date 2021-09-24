if __name__ == '__main__':
    from os.path import join
    # import numpy as np

    from troppo.tasks.core import Task
    from troppo.tasks.task_io import JSONTaskIO

    from cobra.io import read_sbml_model

    base_dir = '/home/scardoso/Documents/PhD/Metabolic_Models'

    '''
    Metabolic Tasks - Essential in All Human Cells for Cell Viability
    '''
    essential_tasks = []
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_1',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorylation of ATP from Glucose',
                                                  'task_group': 'Rephosphorylation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630c': [1, 1000],  # O2
                                                           'MAM01965c': [1, 1000]  # Glucose
                                                           },
                                              outflow_dict={'MAM02040c': [1, 1000],  # H2O
                                                            'MAM01596c': [1, 1000],  # CO2
                                                            'MAM01371s': [1, 1000]
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (0, 1000))},
                                              mandatory_activity=['T_SUCCESS_1_atp_to_adp_task_reaction']
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_2',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorylation of GTP from Glucose',
                                                  'task_group': 'Rephosphorylation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'gtp_to_gdp': ({'MAM02034c': -1, 'MAM02040c': -1,
                                                                             'MAM01948c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_3',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorylation of CTP from Glucose',
                                                  'task_group': 'Rephosphorylation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'ctp_to_cdp': ({'MAM01623c': -1, 'MAM02040c': -1,
                                                                             'MAM01424c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_4',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorylation of UTP from Glucose',
                                                  'task_group': 'Rephosphorylation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'utp_to_udp': ({'MAM03130c': -1, 'MAM02040c': -1,
                                                                             'MAM03106c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_5',
                                              annotations={
                                                  'name': 'ATP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01371c': [1, 1000]  # ATP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_6',
                                              annotations={
                                                  'name': 'CTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01623c': [1, 1000]  # CTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_7',
                                              annotations={
                                                  'name': 'GTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM02034c': [1, 1000]  # GTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_8',
                                              annotations={
                                                  'name': 'UTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM03130c': [1, 1000]  # UTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_9',
                                              annotations={
                                                  'name': 'dATP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01642c': [1, 1000]  # dATP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_10',
                                              annotations={
                                                  'name': 'dCTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01645c': [1, 1000]  # dCTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_11',
                                              annotations={
                                                  'name': 'dGTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01688c': [1, 1000]  # dGTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_12',
                                              annotations={
                                                  'name': 'dTTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01753c': [0, 1000]  # dTTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_13',
                                              annotations={
                                                  'name': 'Histidine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02125s': [1, 1]  # Histidine
                                                           },
                                              outflow_dict={'MAM02125c': [1, 1]  # Histidine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_14',
                                              annotations={
                                                  'name': 'Isoleucine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02184s': [1, 1]  # Isoleucine
                                                           },
                                              outflow_dict={'MAM02184c': [1, 1]  # Isoleucine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_15',
                                              annotations={
                                                  'name': 'Leucine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02360s': [1, 1]  # Leucine
                                                           },
                                              outflow_dict={'MAM02360c': [1, 1]  # Leucine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_16',
                                              annotations={
                                                  'name': 'Lysine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02426s': [1, 1]  # Lysine
                                                           },
                                              outflow_dict={'MAM02426c': [1, 1]  # Lysine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_17',
                                              annotations={
                                                  'name': 'MAMethionine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02471s': [1, 1]  # Methionine
                                                           },
                                              outflow_dict={'MAM02471c': [1, 1]  # Methionine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_18',
                                              annotations={
                                                  'name': 'Phenylalanine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02724s': [1, 1]  # Phenylalanine
                                                           },
                                              outflow_dict={'MAM02724c': [1, 1]  # Phenylalanine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_19',
                                              annotations={
                                                  'name': 'Threonine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02993s': [1, 1]  # Threonine
                                                           },
                                              outflow_dict={'MAM02993c': [1, 1]  # Threonine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_20',
                                              annotations={
                                                  'name': 'Tryptophan Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM03089s': [1, 1]  # Tryptophan
                                                           },
                                              outflow_dict={'MAM03089c': [1, 1]  # Tryptophan
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_21',
                                              annotations={
                                                  'name': 'Valine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM03135s': [1, 1]  # Valine
                                                           },
                                              outflow_dict={'MAM03135c': [1, 1]  # Valine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_22',
                                              annotations={
                                                  'name': 'Glycerate 3-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM00913c': [1, 1000]  # 3-Phospho-D-Glycerate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_23',
                                              annotations={
                                                  'name': 'Mitochondrial acetyl-CoA de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'succinyl_presence': ({'MAM02944m': -1, 'MAM01597m': 1},
                                                                                   (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_24',
                                              annotations={
                                                  'name': 'MAMitochondrial AKG de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01306m': [1, 1000]}  # AKG
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_25',
                                              annotations={
                                                  'name': 'Erythrose 4-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01785c': [1, 1000]  # erythrose 4-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_26',
                                              annotations={
                                                  'name': 'Fructose 6-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01845c': [1, 1000]  # fructose 6-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_27',
                                              annotations={
                                                  'name': 'Glyceraldehyde 3-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01939c': [1, 1000]  # glyceraldehyde 3-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_28',
                                              annotations={
                                                  'name': 'Glucose 6-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01968c': [1, 1000]  # glucose 6-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_29',
                                              annotations={
                                                  'name': 'MAMitochondrial oxaloacetate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM02633m': [1, 1000]  # OAA
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_30',
                                              annotations={
                                                  'name': 'Phosphoenolpyruvate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM02696c': [1, 1000]  # PEP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_31',
                                              annotations={
                                                  'name': 'Pyruvate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM02819c': [1, 1000]  # pyruvate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_32',
                                              annotations={
                                                  'name': 'Ribose 5-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM02845c': [1, 1000]  # ribose 5-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_33',
                                              annotations={
                                                  'name': 'MAMitochondrial succinnyl-CoA de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'succinyl_presence': ({'MAM02944m': -1, 'MAM01597m': 1},
                                                                                   (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_34',
                                              annotations={
                                                  'name': 'Cholesterol de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Other Compounds'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM01450c': [1, 1000]  # cholesterol
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_35',
                                              annotations={
                                                  'name': 'Protein Turnover',
                                                  'task_group': 'Protein Synthesis of Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM01965s': [0, 1000],  # Glucose
                                                           'MAM02578s': [0, 1000],  # NH3
                                                           'MAM02040s': [0, 1000],  # H2O
                                                           'MAM01365s': [0, 1000],  # Arginine
                                                           'MAM02125s': [0, 1000],  # Histidine
                                                           'MAM02426s': [0, 1000],  # Lysine
                                                           'MAM02471s': [0, 1000],  # Methionine
                                                           'MAM02724s': [0, 1000],  # Phenylalanine
                                                           'MAM03089s': [0, 1000],  # Tryptophan
                                                           'MAM03101s': [0, 1000],  # Tyrosine
                                                           'MAM01307s': [0, 1000],  # Alanine
                                                           'MAM01986s': [0, 1000],  # Glycine
                                                           'MAM02896s': [0, 1000],  # Serine
                                                           'MAM02993s': [0, 1000],  # Threonine
                                                           'MAM01370s': [0, 1000],  # Aspartate
                                                           'MAM01974s': [0, 1000],  # Glutamate
                                                           'MAM01369s': [0, 1000],  # Asparagine
                                                           'MAM01975s': [0, 1000],  # Glutamine
                                                           'MAM02184s': [0, 1000],  # Isoleucine
                                                           'MAM02360s': [0, 1000],  # Leucine
                                                           'MAM02770s': [0, 1000],  # Proline
                                                           'MAM03135s': [0, 1000],  # Valine
                                                           'MAM01628s': [0, 1000]  # Cysteine
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM02042s': [0, 1000],  # H2S
                                                            'MAM03121s': [0, 1000],  # urea
                                                            'MAM01308c': [0.0001, 1000]  # albumin
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_36',
                                              annotations={
                                                  'name': 'Oxidative Phosphorylation',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM02553m': [1, 1],  # NADH
                                                           'MAM02943m': [1, 1]  # succinate
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM02552m': [1, 1],  # NAD+
                                                            'MAM01862m': [1, 1]  # fumarate
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_37',
                                              annotations={
                                                  'name': 'Oxidative Decarboxylation',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'MAM02819m': [1, 1],  # pyruvate
                                                           'MAM02552m': [1, 1],  # NAD+
                                                           'MAM01597m': [1, 1]  # CoA
                                                           },
                                              outflow_dict={'MAM01261m': [1, 1],  # acetyl-CoA
                                                            'MAM02553m': [1, 1],  # NADH
                                                            'MAM01596s': [1, 1],  # CO2
                                                            'MAM02039m': [0, 1000]  # H+
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_38',
                                              annotations={
                                                  'name': 'Krebs Cycle NADH',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'MAM01261m': [1, 1],  # acetyl-CoA
                                                           'MAM01948m': [1, 1],  # GDP
                                                           'MAM03103m': [1, 1],  # ubiquinone
                                                           'MAM02552m': [3, 3],  # NAD+
                                                           'MAM02751m': [0, 1000],  # Pi
                                                           'MAM02040s': [0, 1000]  # H2O
                                                           },
                                              outflow_dict={'MAM01597m': [1, 1],  # CoA
                                                            'MAM03102m': [1, 1],  # ubiquinol
                                                            'MAM02034m': [1, 1],  # GTP
                                                            'MAM02553m': [3, 3],  # NADH
                                                            'MAM01596s': [0, 1000],  # CO2
                                                            'MAM02039c': [0, 1000],  # H+
                                                            'MAM02039m': [0, 1000]  # H+
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_39',
                                              annotations={
                                                  'name': 'Ubiquinol-to-proton',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'MAM03102m': [1, 1],  # ubiquinol
                                                           'MAM02039m': [0, 1000],  # H+
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM03103m': [1, 1],  # ubiquinone
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM02039c': [6, 1000]  # H+
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_40',
                                              annotations={
                                                  'name': 'Ubiquinol-to-ATP',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'MAM03102m': [1, 1],  # ubiquinol
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM03103m': [1, 1],  # ubiquinone
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM02039c': [0, 1000],  # H+
                                                            'MAM02039m': [0, 1000]  # H+
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_41',
                                              annotations={
                                                  'name': 'Beta Oxidation of Saturated FA',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02938s': [1, 1],  # stearate
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_42',
                                              annotations={
                                                  'name': 'Beta Oxidation of Long-chain FA',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM00315s': [1, 1],  # 12,15,18,21-tetracosatetraenoic acid
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_43',
                                              annotations={
                                                  'name': 'Beta Oxidation of Odd-chain FA',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02456s': [1, 1],  # margaric acid
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_44',
                                              annotations={
                                                  'name': 'Beta Oxidation of Unsaturated FA (n-9)',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM03153s': [1, 1],  # ximenic acid
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_45',
                                              annotations={
                                                  'name': 'Beta Oxidation of Unsaturated FA (n-6)',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02387s': [1, 1],  # Linoleate
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_46',
                                              annotations={
                                                  'name': 'Uptake and Beta Oxidation of all NEFAs',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02560s': [1, 1],  # NEFA blood pool in
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (0, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_47',
                                              annotations={
                                                  'name': 'Choline Uptake',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'MAM01513s': [1, 1]  # Choline
                                                           },
                                              outflow_dict={'MAM01513c': [1, 1]  # Choline
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_48',
                                              annotations={
                                                  'name': 'Inositol Uptake',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02171s': [1, 1]  # Inositol
                                                           },
                                              outflow_dict={'MAM02171c': [1, 1]  # Inositol
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_49',
                                              annotations={
                                                  'name': 'Phosphatidylcholine de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'MAM01513s': [0, 1000],  # choline
                                                           'MAM01965s': [0, 1000],  # glucose
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM02560s': [0, 1000],  # NEFA blood pool in
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02684c': [1, 1000],  # PC-LD pool
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (-1000, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_50',
                                              annotations={
                                                  'name': 'Phosphatidylethanolamine de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'MAM01797s': [0, 1000],  # ethanolamine
                                                           'MAM01965s': [0, 1000],  # glucose
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM02560s': [0, 1000],  # NEFA blood pool in
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02685c': [1, 1000],  # PE-LD pool
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_51',
                                              annotations={
                                                  'name': 'Phosphatidylserine de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02896s': [0, 1000],  # serine
                                                           'MAM01965s': [0, 1000],  # glucose
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM02560s': [0, 1000],  # NEFA blood pool in
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02808c': [1, 1000],  # PS-LD pool
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_52',
                                              annotations={
                                                  'name': 'Phosphatidylinositol de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'MAM02171s': [0, 1000],  # inositol
                                                           'MAM01965s': [0, 1000],  # glucose
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM02560s': [0, 1000],  # NEFA blood pool in
                                                           'MAM02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'MAM02750c': [1, 1000],  # Pi pool
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_53',
                                              annotations={
                                                  'name': 'Thiamin Phosphorylation to TPP',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'MAM02982s': [0, 1000],  # thiamin
                                                           'MAM02751s': [0, 1000],  # Pi
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM02040s': [0, 1000],  # H2O
                                                           'MAM02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02984c': [1, 1000],  # thiamin-PP
                                                            'MAM02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (-1000, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_54',
                                              annotations={
                                                  'name': 'Coenzyme A Synthesis from Pantothenate',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'MAM02680s': [0, 1000],  # Pantothenate
                                                           'MAM01628s': [0, 1000],  # cysteine
                                                           'MAM01965s': [0, 1000],  # glucose
                                                           'MAM02751s': [0, 1000],  # Pi
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM02040s': [0, 1000],  # H2O
                                                           'MAM02578s': [0, 1000]  # NH3
                                                           },
                                              outflow_dict={'MAM01597c': [1, 1000],  # CoA
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1, 'MAM02039c': 1},
                                                                            (-1000, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_55',
                                              annotations={
                                                  'name': 'FAD Synthesis from Riboflavin',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'MAM01965s': [0, 1000],  # glucose
                                                           'MAM02751s': [0, 1000],  # Pi
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM02842s': [0, 1000],  # riboflavin
                                                           'MAM02578s': [0, 1000]  # NH3
                                                           },
                                              outflow_dict={'MAM01802s': [1, 1000],  # FAD
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_56',
                                              annotations={
                                                  'name': 'Heme Biosynthesis',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'MAM01965s': [0, 1000],  # glucose
                                                           'MAM02630s': [0, 1000],  # O2
                                                           'MAM01821s': [0, 1000],  # Fe2+
                                                           'MAM02578s': [0, 1000]  # NH3
                                                           },
                                              outflow_dict={'MAM02049c': [1, 1000],  # heme
                                                            'MAM02040s': [0, 1000],  # H2O
                                                            'MAM01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_1',
                                              annotations={
                                                  'name': 'Oxygen from Water',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02040s': [0, 1000]  # H2O
                                                           },
                                              outflow_dict={'MAM02630s': [1, 1000],  # O2
                                                            'MAM02039s': [0, 1000]  # H+
                                                            },
                                              reaction_dict={'h_transport_cs': ({'MAM02039c': -1, 'MAM02039s': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_mc': ({'MAM02039m': -1, 'MAM02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_pc': ({'MAM02039p': -1, 'MAM02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_gc': ({'MAM02039g': -1, 'MAM02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_rc': ({'MAM02039r': -1, 'MAM02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_lc': ({'MAM02039l': -1, 'MAM02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_nc': ({'MAM02039n': -1, 'MAM02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_ic': ({'MAM02039i': -1, 'MAM02039c': 1},
                                                                                (-1000, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_2',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Cytoplasm',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02040c': [0, 1000],  # H2O
                                                           'MAM02630c': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040c': [0, 1000],  # H2O
                                                            'MAM02630c': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1,
                                                                             'MAM02039c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_3',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Mitochondria',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02040m': [0, 1000],  # H2O
                                                           'MAM02630m': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040m': [0, 1000],  # H2O
                                                            'MAM02630m': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371m': -1, 'MAM02040m': -1,
                                                                             'MAM01285m': 1, 'MAM02751m': 1,
                                                                             'MAM02039m': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_4',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Peroxisome',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02040p': [0, 1000],  # H2O
                                                           'MAM02630p': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040p': [0, 1000],  # H2O
                                                            'MAM02630p': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371p': -1, 'MAM02040p': -1,
                                                                             'MAM01285p': 1, 'MAM02751p': 1,
                                                                             'MAM02039p': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    #essential_tasks = essential_tasks + [Task(name='T_FAIL_5',
    #                                          annotations={
    #                                              'name': 'Rephosphorylation of ATP in Golgi',
    #                                              'task_group': 'Should Fail'},
    #                                          should_fail=True,
    #                                          inflow_dict={'MAM02040g': [0, 1000],  # H2O
    #                                                       'MAM02630g': [0, 1000]  # O2
    #                                                       },
    #                                          outflow_dict={'MAM02040g': [0, 1000],  # H2O
    #                                                        'MAM02630g': [0, 1000]  # O2
    #                                                        },
    #                                          reaction_dict={'atp_to_adp': ({'MAM01371g': -1, 'MAM02040g': -1,
    #                                                                         'MAM01285g': 1, 'MAM02751g': 1,
    #                                                                         'MAM02039g': 1},
    #                                                                        (1, 1000))
    #                                                         }
    #                                          )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_6',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Endoplasmatic Reticulum',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02040r': [0, 1000],  # H2O
                                                           'MAM02630r': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040r': [0, 1000],  # H2O
                                                            'MAM02630r': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371r': -1, 'MAM02040r': -1,
                                                                             'MAM01285r': 1, 'MAM02751r': 1,
                                                                             'MAM02039r': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_7',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Nucleus',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02040n': [0, 1000],  # H2O
                                                           'MAM02630n': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'MAM02040n': [0, 1000],  # H2O
                                                            'MAM02630n': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371n': -1, 'MAM02040n': -1,
                                                                             'MAM01285n': 1, 'MAM02751n': 1,
                                                                             'MAM02039n': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    # essential_tasks = essential_tasks + [Task(name='T_FAIL_8',
    #                                          annotations={
    #                                              'name': 'Rephosphorylation of ATP in Lysosome',
    #                                              'task_group': 'Should Fail'},
    #                                          should_fail=True,
    #                                          inflow_dict={'MAM02040l': [0, 1000],  # H2O
    #                                                       'MAM02630l': [0, 1000]  # O2
    #                                                       },
    #                                          outflow_dict={'MAM02040l': [0, 1000],  # H2O
    #                                                        'MAM02630l': [0, 1000]  # O2
    #                                                        },
    #                                          reaction_dict={'atp_to_adp': ({'MAM01371l': -1, 'MAM02040l': -1,
    #                                                                         'MAM01285l': 1, 'MAM02751l': 1,
    #                                                                         'MAM02039l': 1},
    #                                                                        (1, 1000))
    #                                                         }
    #                                          )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_9',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP from Protons',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02040c': [0, 1000],  # H2O
                                                           'MAM02630c': [0, 1000],  # O2
                                                           'MAM02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02040c': [0, 1000],  # H2O
                                                            'MAM02630c': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1,
                                                                             'MAM02039c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_10',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Cytoplasm',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1,
                                                                             'MAM02039c': 1},
                                                                            (-1000, 1000)),
                                                             'nadh_to_nad': ({'MAM02553c': -2, 'MAM02630c': -1,
                                                                              'MAM02039c': -2,
                                                                              'MAM02552c': 2, 'MAM02040c': 2},
                                                                             (1, 1000)),
                                                             'nadph_to_nadp': ({'MAM02555c': -1, 'MAM02552c': -1,
                                                                                'MAM02554c': 1, 'MAM02553c': 1},
                                                                               (-1000, 1000)),
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_11',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Mitochondria',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM02039m': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'atp_to_adp': ({'MAM01371m': -1, 'MAM02040m': -1,
                                                                             'MAM01285m': 1, 'MAM02751m': 1,
                                                                             'MAM02039m': 1},
                                                                            (-1000, 1000)),
                                                             'nadh_to_nad': ({'MAM02553m': -2, 'MAM02630m': -1,
                                                                              'MAM02039m': -2,
                                                                              'MAM02552m': 2, 'MAM02040m': 2},
                                                                             (1, 1000)),
                                                             'nadph_to_nadp': ({'MAM02555m': -1, 'MAM02552m': -1,
                                                                                'MAM02554m': 1, 'MAM02553m': 1},
                                                                               (-1000, 1000)),
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_12',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Peroxisome',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM02039p': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'nadh_to_nad': ({'MAM02553m': -2, 'MAM02630m': -1,
                                                                              'MAM02039m': -2,
                                                                              'MAM02552m': 2, 'MAM02040m': 2},
                                                                             (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_13',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Cytoplasm (FADH2)',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'fadh2_to_fad': ({'MAM01803c': -1, 'MAM02630c': -1,
                                                                               'MAM01802c': 1, 'MAM02040c': 1},
                                                                              (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_14',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Mitochondria (FADH2)',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM02039m': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'fadh2_to_fad': ({'MAM01803m': -1, 'MAM02630m': -1,
                                                                               'MAM01802m': 1, 'MAM02040m': 1},
                                                                              (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_15',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Peroxisome (FADH2)',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02630s': [0, 1000],  # O2
                                                           'MAM02039p': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'MAM02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'fadh2_to_fad': ({'MAM01803p': -1, 'MAM02630p': -1,
                                                                               'MAM01802p': 1, 'MAM02040p': 1},
                                                                              (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_16',
                                              annotations={
                                                  'name': 'Generation of CO2 from nothing',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={},
                                              outflow_dict={'MAM01596s': [0.001, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_17',
                                              annotations={
                                                  'name': 'ATP  Phosphorilated per mol of Glucose Exceeds 40 mol',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM01965s': [0, 1],  # Glucose
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={},
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1,
                                                                             'MAM02039c': 1},
                                                                            (40, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_18',
                                              annotations={
                                                  'name': 'Free Transport of phosphate from cytoplasm'
                                                          'to mitochondria',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={},
                                              outflow_dict={},
                                              reaction_dict={'pim_to_pic': ({'MAM02751m': -1, 'MAM02751c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_19',
                                              annotations={
                                                  'name': 'Anaerobic Production of Propanoate from Glucose',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM01965s': [0, 1]},  # Glucose
                                              outflow_dict={'MAM02772s': [1, 1000]},  # Propanoate
                                              reaction_dict={}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_20',
                                              annotations={
                                                  'name': 'Anaerobic ATP phosphorylation per glucose exceeds 2',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM01965s': [0, 1]},  # Glucose
                                              outflow_dict={},
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1,
                                                                             'MAM02039c': 1},
                                                                            (2.5, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_21',
                                              annotations={
                                                  'name': 'ATP Production from Pi and O2',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'MAM02751s': [0, 1000],  # Pi
                                                           'MAM02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={},
                                              reaction_dict={'atp_to_adp': ({'MAM01371c': -1, 'MAM02040c': -1,
                                                                             'MAM01285c': 1, 'MAM02751c': 1,
                                                                             'MAM02039c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    jtio = JSONTaskIO()
    jtio.write_task(join(base_dir, 'GENERAL/utility_data/metabolic_tasks_hsa_cellViability_consistent.json'),
                    essential_tasks)

    # --- Metabolic Tasks for all Human Cells - for HumanGEM_consistent
    # The following tasks were removed:
    # T_FAIL_5, T_FAIL_8

    '''
    HumanGEM GPRs and GENES
    '''

    # Directories:
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    general_utility_data_dir = join(base_dir, 'GENERAL/utility_data')

    # --- Read HumanGEM-1.8.0 original model ---
    print('\nReading HumanGEM-1.8.0...')
    HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.8.0.xml.gz'))
    # --- Get HumanGEM model's GPRs ---
    print('Collecting GPR rules from HumanGEM-1.8.0')
    reaction_ids = []
    for i in HumanGEM.reactions:
        reaction_ids.append(i.id)
    f = open(join(general_utility_data_dir, 'HumanGEM-1.8.0_GPRs.txt'), 'w')
    for item in reaction_ids:
        f.write("%s\t%s\n" % (item, HumanGEM.reactions.get_by_id(item).gene_reaction_rule))
    f.close()
    # --- Get HumanGEM model's genes ---
    print('Collecting genes from HumanGEM-1.8.0')
    gene_ids = []
    for i in HumanGEM.genes:
        gene_ids.append(i.id)
    f = open(join(general_utility_data_dir, 'HumanGEM-1.8.0_GENES.txt'), 'w')
    for item in gene_ids:
        f.write("%s\n" % item)
    f.close()

    # --- Read HumanGEM_1.4.1_consistent original model ---
    print('\nReading HumanGEM-1.8.0_consistent...')
    HumanGEM_consistent = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.8.0_consistent.xml.gz'))
    # --- Get HumanGEM model's GPRs ---
    print('Collecting GPR rules from HumanGEM-1.8.0_consistent')
    reaction_ids = []
    for i in HumanGEM_consistent.reactions:
        reaction_ids.append(i.id)
    f = open(join(general_utility_data_dir, 'HumanGEM-1.8.0_consistent_GPRs.txt'), 'w')
    for item in reaction_ids:
        f.write("%s\t%s\n" % (item, HumanGEM_consistent.reactions.get_by_id(item).gene_reaction_rule))
    f.close()
    # --- Get HumanGEM model's genes ---
    print('Collecting genes from HumanGEM-1.8.0_consistent')
    gene_ids = []
    for i in HumanGEM_consistent.genes:
        gene_ids.append(i.id)
    f = open(join(general_utility_data_dir, 'HumanGEM-1.8.0_consistent_GENES.txt'), 'w')
    for item in gene_ids:
        f.write("%s\n" % item)
    f.close()
