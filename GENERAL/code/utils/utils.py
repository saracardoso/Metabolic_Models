if __name__ == '__main__':
    from os.path import join
    import numpy as np

    from troppo.tasks.core import Task
    from troppo.tasks.task_io import JSONTaskIO

    base_dir = '/home/scardoso/Documents/PhD'

    '''
    Metabolic Tasks - Essential in All Human Cells for Cell Viability
    '''
    essential_tasks = []
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_1',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorilation of ATP from Glucose',
                                                  'task_group': 'Rephosphorilation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_2',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorilation of GTP from Glucose',
                                                  'task_group': 'Rephosphorilation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'gtp_to_gdp': ({'m02034c': -1, 'm02040c': -1,
                                                                             'm01948c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_3',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorilation of CTP from Glucose',
                                                  'task_group': 'Rephosphorilation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'ctp_to_cdp': ({'m01623c': -1, 'm02040c': -1,
                                                                             'm01424c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_4',
                                              annotations={
                                                  'name': 'Aerobic Rephosphorilation of UTP from Glucose',
                                                  'task_group': 'Rephosphorilation of Nucleoside Triphosphates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'utp_to_udp': ({'m03130c': -1, 'm02040c': -1,
                                                                             'm03106c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_5',
                                              annotations={
                                                  'name': 'ATP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01371c': [1, 1000]  # ATP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_6',
                                              annotations={
                                                  'name': 'CTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01623c': [1, 1000]  # CTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_7',
                                              annotations={
                                                  'name': 'GTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm02034c': [1, 1000]  # GTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_8',
                                              annotations={
                                                  'name': 'UTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm03130c': [1, 1000]  # UTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_9',
                                              annotations={
                                                  'name': 'dATP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01642c': [1, 1000]  # dATP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_10',
                                              annotations={
                                                  'name': 'dCTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01645c': [1, 1000]  # dCTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_11',
                                              annotations={
                                                  'name': 'dGTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01688c': [1, 1000]  # dGTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_12',
                                              annotations={
                                                  'name': 'dTTP de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Nucleotides'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01753c': [0, 1000]  # dTTP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_13,',
                                              annotations={
                                                  'name': 'Histidine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02125s': [1, 1]  # Histidine
                                                           },
                                              outflow_dict={'m02125c': [1, 1]  # Histidine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_14',
                                              annotations={
                                                  'name': 'Isoleucine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02184s': [1, 1]  # Isoleucine
                                                           },
                                              outflow_dict={'m02184c': [1, 1]  # Isoleucine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_15',
                                              annotations={
                                                  'name': 'Leucine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02360s': [1, 1]  # Leucine
                                                           },
                                              outflow_dict={'m02360c': [1, 1]  # Leucine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_16',
                                              annotations={
                                                  'name': 'Lysine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02426s': [1, 1]  # Lysine
                                                           },
                                              outflow_dict={'m02426c': [1, 1]  # Lysine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_17',
                                              annotations={
                                                  'name': 'Methionine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02471s': [1, 1]  # Methionine
                                                           },
                                              outflow_dict={'m02471c': [1, 1]  # Methionine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_18',
                                              annotations={
                                                  'name': 'Phenylalanine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02724s': [1, 1]  # Phenylalanine
                                                           },
                                              outflow_dict={'m02724c': [1, 1]  # Phenylalanine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_19',
                                              annotations={
                                                  'name': 'Threonine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02993s': [1, 1]  # Threonine
                                                           },
                                              outflow_dict={'m02993c': [1, 1]  # Threonine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_20',
                                              annotations={
                                                  'name': 'Tryptophan Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m03089s': [1, 1]  # Tryptophan
                                                           },
                                              outflow_dict={'m03089c': [1, 1]  # Tryptophan
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_21',
                                              annotations={
                                                  'name': 'Valine Uptake',
                                                  'task_group': 'Uptake of Essential Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m03135s': [1, 1]  # Valine
                                                           },
                                              outflow_dict={'m03135c': [1, 1]  # Valine
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_22',
                                              annotations={
                                                  'name': 'Glycerate 3-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm00913c': [1, 1000]  # 3-Phospho-D-Glycerate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_23',
                                              annotations={
                                                  'name': 'Mitochondrial acetyl-CoA de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'succinyl_presence': ({'m02944m': -1, 'm01597m': 1},
                                                                                   (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_24',
                                              annotations={
                                                  'name': 'Mitochondrial AKG de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01306m': [1, 1000]}  # AKG
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_25',
                                              annotations={
                                                  'name': 'Erythrose 4-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01785c': [1, 1000]  # erythrose 4-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_26',
                                              annotations={
                                                  'name': 'Fructose 6-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01845c': [1, 1000]  # fructose 6-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_27',
                                              annotations={
                                                  'name': 'Glyceraldehyde 3-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01939c': [1, 1000]  # glyceraldehyde 3-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_28',
                                              annotations={
                                                  'name': 'Glucose 6-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01968c': [1, 1000]  # glucose 6-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_29',
                                              annotations={
                                                  'name': 'Mitochondrial oxaloacetate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm02633m': [1, 1000]  # OAA
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_30',
                                              annotations={
                                                  'name': 'Phosphoenolpyruvate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm02696c': [1, 1000]  # PEP
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_31',
                                              annotations={
                                                  'name': 'Pyruvate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm02819c': [1, 1000]  # pyruvate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_32',
                                              annotations={
                                                  'name': 'Ribose 5-phosphate de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm02845c': [1, 1000]  # ribose 5-phosphate
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_33',
                                              annotations={
                                                  'name': 'Mitochondrial succinnyl-CoA de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Key Intermediates'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'succinyl_presence': ({'m02944m': -1, 'm01597m': 1},
                                                                                   (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_34',
                                              annotations={
                                                  'name': 'Cholesterol de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Other Compounds'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000]  # Glucose
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm01450c': [1, 1000]  # cholesterol
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_35',
                                              annotations={
                                                  'name': 'Protein Turnover',
                                                  'task_group': 'Protein Synthesis of Amino acids'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm01965s': [0, 1000],  # Glucose
                                                           'm02578s': [0, 1000],  # NH3
                                                           'm02040s': [0, 1000],  # H2O
                                                           'm01365s': [0, 1000],  # Arginine
                                                           'm02125s': [0, 1000],  # Histidine
                                                           'm02426s': [0, 1000],  # Lysine
                                                           'm02471s': [0, 1000],  # Methionine
                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                           'm03089s': [0, 1000],  # Tryptophan
                                                           'm03101s': [0, 1000],  # Tyrosine
                                                           'm01307s': [0, 1000],  # Alanine
                                                           'm01986s': [0, 1000],  # Glycine
                                                           'm02896s': [0, 1000],  # Serine
                                                           'm02993s': [0, 1000],  # Threonine
                                                           'm01370s': [0, 1000],  # Aspartate
                                                           'm01974s': [0, 1000],  # Glutamate
                                                           'm01369s': [0, 1000],  # Asparagine
                                                           'm01975s': [0, 1000],  # Glutamine
                                                           'm02184s': [0, 1000],  # Isoleucine
                                                           'm02360s': [0, 1000],  # Leucine
                                                           'm02770s': [0, 1000],  # Proline
                                                           'm03135s': [0, 1000],  # Valine
                                                           'm01628s': [0, 1000]  # Cysteine
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm02042s': [0, 1000],  # H2S
                                                            'm03121s': [0, 1000],  # urea
                                                            'm01308c': [0.0001, 1000]  # albumin
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_36',
                                              annotations={
                                                  'name': 'Oxidative Phosphorylation',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm02553m': [1, 1],  # NADH
                                                           'm02943m': [1, 1]  # succinate
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm02552m': [1, 1],  # NAD+
                                                            'm01862m': [1, 1]  # fumarate
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_37',
                                              annotations={
                                                  'name': 'Oxidative Decarboxylation',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'m02819m': [1, 1],  # pyruvate
                                                           'm02552m': [1, 1],  # NAD+
                                                           'm01597m': [1, 1]  # CoA
                                                           },
                                              outflow_dict={'m01261m': [1, 1],  # acetyl-CoA
                                                            'm02553m': [1, 1],  # NADH
                                                            'm01596s': [1, 1],  # CO2
                                                            'm02039m': [0, 1000]  # H+
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_38',
                                              annotations={
                                                  'name': 'Krebs Cycle NADH',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'m01261m': [1, 1],  # acetyl-CoA
                                                           'm01948m': [1, 1],  # GDP
                                                           'm03103m': [1, 1],  # ubiquinone
                                                           'm02552m': [3, 3],  # NAD+
                                                           'm02751m': [0, 1000],  # Pi
                                                           'm02040s': [0, 1000]  # H2O
                                                           },
                                              outflow_dict={'m01597m': [1, 1],  # CoA
                                                            'm03102m': [1, 1],  # ubiquinol
                                                            'm02034m': [1, 1],  # GTP
                                                            'm02553m': [3, 3],  # NADH
                                                            'm01596s': [0, 1000],  # CO2
                                                            'm02039c': [0, 1000],  # H+
                                                            'm02039m': [0, 1000]  # H+
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_39',
                                              annotations={
                                                  'name': 'Ubiquinol-to-proton',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'m03102m': [1, 1],  # ubiquinol
                                                           'm02039m': [0, 1000],  # H+
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m03103m': [1, 1],  # ubiquinone
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm02039c': [6, 1000]  # H+
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_40',
                                              annotations={
                                                  'name': 'Ubiquinol-to-ATP',
                                                  'task_group': 'Electron Transport Chain and TCA'},
                                              should_fail=False,
                                              inflow_dict={'m03102m': [1, 1],  # ubiquinol
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m03103m': [1, 1],  # ubiquinone
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm02039c': [0, 1000],  # H+
                                                            'm02039m': [0, 1000]  # H+
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (1, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_41',
                                              annotations={
                                                  'name': 'Beta Oxidation of Saturated FA',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'m02938s': [1, 1],  # stearate
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_42',
                                              annotations={
                                                  'name': 'Beta Oxidation of Long-chain FA',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'m00315s': [1, 1],  # 12,15,18,21-tetracosatetraenoic acid
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_43',
                                              annotations={
                                                  'name': 'Beta Oxidation of Odd-chain FA',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'m02456s': [1, 1],  # margaric acid
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_44',
                                              annotations={
                                                  'name': 'Beta Oxidation of Unsaturated FA (n-9)',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'m03153s': [1, 1],  # ximenic acid
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_45',
                                              annotations={
                                                  'name': 'Beta Oxidation of Unsaturated FA (n-6)',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'m02387s': [1, 1],  # Linoleate
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_46',
                                              annotations={
                                                  'name': 'Uptake and Beta Oxidation of all NEFAs',
                                                  'task_group': 'Beta Oxidation of Fatty Acids'},
                                              should_fail=False,
                                              inflow_dict={'m02560s': [1, 1],  # NEFA blood pool in
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (0, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_47',
                                              annotations={
                                                  'name': 'Choline Uptake',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'m01513s': [1, 1]  # Choline
                                                           },
                                              outflow_dict={'m01513c': [1, 1]  # Choline
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_48',
                                              annotations={
                                                  'name': 'Inositol Uptake',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'m02171s': [1, 1]  # Inositol
                                                           },
                                              outflow_dict={'m02171c': [1, 1]  # Inositol
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_49',
                                              annotations={
                                                  'name': 'Phosphatidylcholine de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'m01513s': [0, 1000],  # choline
                                                           'm01965s': [0, 1000],  # glucose
                                                           'm02630s': [0, 1000],  # O2
                                                           'm02560s': [0, 1000],  # NEFA blood pool in
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02684c': [1, 1000],  # PC-LD pool
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (-1000, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_50',
                                              annotations={
                                                  'name': 'Phosphatidylethanolamine de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'m01797s': [0, 1000],  # ethanolamine
                                                           'm01965s': [0, 1000],  # glucose
                                                           'm02630s': [0, 1000],  # O2
                                                           'm02560s': [0, 1000],  # NEFA blood pool in
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02685c': [1, 1000],  # PE-LD pool
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_51',
                                              annotations={
                                                  'name': 'Phosphatidylserine de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'m02896s': [0, 1000],  # serine
                                                           'm01965s': [0, 1000],  # glucose
                                                           'm02630s': [0, 1000],  # O2
                                                           'm02560s': [0, 1000],  # NEFA blood pool in
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02808c': [1, 1000],  # PS-LD pool
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_52',
                                              annotations={
                                                  'name': 'Phosphatidylinositol de novo Synthesis',
                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                              should_fail=False,
                                              inflow_dict={'m02171s': [0, 1000],  # inositol
                                                           'm01965s': [0, 1000],  # glucose
                                                           'm02630s': [0, 1000],  # O2
                                                           'm02560s': [0, 1000],  # NEFA blood pool in
                                                           'm02751s': [0, 1000]  # Pi
                                                           },
                                              outflow_dict={'m02750c': [1, 1000],  # Pi pool
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_53',
                                              annotations={
                                                  'name': 'Thiamin Phosphorylation to TPP',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'m02982s': [0, 1000],  # thiamin
                                                           'm02751s': [0, 1000],  # Pi
                                                           'm02630s': [0, 1000],  # O2
                                                           'm02040s': [0, 1000],  # H2O
                                                           'm02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02984c': [1, 1000],  # thiamin-PP
                                                            'm02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (-1000, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_54',
                                              annotations={
                                                  'name': 'Coenzyme A Synthesis from Pantothenate',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'m02680s': [0, 1000],  # Pantothenate
                                                           'm01628s': [0, 1000],  # cysteine
                                                           'm01965s': [0, 1000],  # glucose
                                                           'm02751s': [0, 1000],  # Pi
                                                           'm02630s': [0, 1000],  # O2
                                                           'm02040s': [0, 1000],  # H2O
                                                           'm02578s': [0, 1000]  # NH3
                                                           },
                                              outflow_dict={'m01597c': [1, 1000],  # CoA
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1, 'm02039c': 1},
                                                                            (-1000, 1000))}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_55',
                                              annotations={
                                                  'name': 'FAD Synthesis from Riboflavin',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'m01965s': [0, 1000],  # glucose
                                                           'm02751s': [0, 1000],  # Pi
                                                           'm02630s': [0, 1000],  # O2
                                                           'm02842s': [0, 1000],  # riboflavin
                                                           'm02578s': [0, 1000]  # NH3
                                                           },
                                              outflow_dict={'m01802s': [1, 1000],  # FAD
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_56',
                                              annotations={
                                                  'name': 'Heme Biosynthesis',
                                                  'task_group': 'Vitamins and co-factors'},
                                              should_fail=False,
                                              inflow_dict={'m01965s': [0, 1000],  # glucose
                                                           'm02630s': [0, 1000],  # O2
                                                           'm01821s': [0, 1000],  # Fe2+
                                                           'm02578s': [0, 1000]  # NH3
                                                           },
                                              outflow_dict={'m02049c': [1, 1000],  # heme
                                                            'm02040s': [0, 1000],  # H2O
                                                            'm01596s': [0, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_1',
                                              annotations={
                                                  'name': 'Oxygen from Water',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040s': [0, 1000]  # H2O
                                                           },
                                              outflow_dict={'m02630s': [1, 1000],  # O2
                                                            'm02039s': [0, 1000]  # H+
                                                            },
                                              reaction_dict={'h_transport_cs': ({'m02039c': -1, 'm02039s': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_mc': ({'m02039m': -1, 'm02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_pc': ({'m02039p': -1, 'm02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_gc': ({'m02039g': -1, 'm02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_rc': ({'m02039r': -1, 'm02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_lc': ({'m02039l': -1, 'm02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_nc': ({'m02039n': -1, 'm02039c': 1},
                                                                                (-1000, 1000)),
                                                             'h_transport_ic': ({'m02039i': -1, 'm02039c': 1},
                                                                                (-1000, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_2',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Cytoplasm',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040c': [0, 1000],  # H2O
                                                           'm02630c': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040c': [0, 1000],  # H2O
                                                            'm02630c': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1,
                                                                             'm02039c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_3',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Mitochondria',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040m': [0, 1000],  # H2O
                                                           'm02630m': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040m': [0, 1000],  # H2O
                                                            'm02630m': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371m': -1, 'm02040m': -1,
                                                                             'm01285m': 1, 'm02751m': 1,
                                                                             'm02039m': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_4',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Peroxisome',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040p': [0, 1000],  # H2O
                                                           'm02630p': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040p': [0, 1000],  # H2O
                                                            'm02630p': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371p': -1, 'm02040p': -1,
                                                                             'm01285p': 1, 'm02751p': 1,
                                                                             'm02039p': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_5',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Golgi',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040g': [0, 1000],  # H2O
                                                           'm02630g': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040g': [0, 1000],  # H2O
                                                            'm02630g': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371g': -1, 'm02040g': -1,
                                                                             'm01285g': 1, 'm02751g': 1,
                                                                             'm02039g': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_6',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Endoplasmatic Reticulum',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040r': [0, 1000],  # H2O
                                                           'm02630r': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040r': [0, 1000],  # H2O
                                                            'm02630r': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371r': -1, 'm02040r': -1,
                                                                             'm01285r': 1, 'm02751r': 1,
                                                                             'm02039r': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_7',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Nucleus',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040n': [0, 1000],  # H2O
                                                           'm02630n': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040n': [0, 1000],  # H2O
                                                            'm02630n': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371n': -1, 'm02040n': -1,
                                                                             'm01285n': 1, 'm02751n': 1,
                                                                             'm02039n': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_8',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP in Lysosome',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040l': [0, 1000],  # H2O
                                                           'm02630l': [0, 1000]  # O2
                                                           },
                                              outflow_dict={'m02040l': [0, 1000],  # H2O
                                                            'm02630l': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371l': -1, 'm02040l': -1,
                                                                             'm01285l': 1, 'm02751l': 1,
                                                                             'm02039l': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_9',
                                              annotations={
                                                  'name': 'Rephosphorylation of ATP from Protons',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02040c': [0, 1000],  # H2O
                                                           'm02630c': [0, 1000],  # O2
                                                           'm02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02040c': [0, 1000],  # H2O
                                                            'm02630c': [0, 1000]  # O2
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1,
                                                                             'm02039c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_10',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Cytoplasm',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1,
                                                                             'm02039c': 1},
                                                                            (-1000, 1000)),
                                                             'nadh_to_nad': ({'m02553c': -2, 'm02630c': -1,
                                                                              'm02039c': -2,
                                                                              'm02552c': 2, 'm02040c': 2},
                                                                             (1, 1000)),
                                                             'nadph_to_nadp': ({'m02555c': -1, 'm02552c': -1,
                                                                                'm02554c': 1, 'm02553c': 1},
                                                                               (-1000, 1000)),
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_11',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Mitochondria',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm02039m': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'atp_to_adp': ({'m01371m': -1, 'm02040m': -1,
                                                                             'm01285m': 1, 'm02751m': 1,
                                                                             'm02039m': 1},
                                                                            (-1000, 1000)),
                                                             'nadh_to_nad': ({'m02553m': -2, 'm02630m': -1,
                                                                              'm02039m': -2,
                                                                              'm02552m': 2, 'm02040m': 2},
                                                                             (1, 1000)),
                                                             'nadph_to_nadp': ({'m02555m': -1, 'm02552m': -1,
                                                                                'm02554m': 1, 'm02553m': 1},
                                                                               (-1000, 1000)),
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_12',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Peroxisome',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm02039p': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'nadh_to_nad': ({'m02553m': -2, 'm02630m': -1,
                                                                              'm02039m': -2,
                                                                              'm02552m': 2, 'm02040m': 2},
                                                                             (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_13',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Cytoplasm (FADH2)',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm02039c': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'fadh2_to_fad': ({'m01803c': -1, 'm02630c': -1,
                                                                               'm01802c': 1, 'm02040c': 1},
                                                                              (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_14',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Mitochondria (FADH2)',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm02039m': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'fadh2_to_fad': ({'m01803m': -1, 'm02630m': -1,
                                                                               'm01802m': 1, 'm02040m': 1},
                                                                              (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_15',
                                              annotations={
                                                  'name': 'Generation of Reduction Potential in Peroxisome (FADH2)',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                           'm02039p': [0, 1000]  # H+
                                                           },
                                              outflow_dict={'m02040s': [0, 1000]  # H2O
                                                            },
                                              reaction_dict={'fadh2_to_fad': ({'m01803p': -1, 'm02630p': -1,
                                                                               'm01802p': 1, 'm02040p': 1},
                                                                              (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_16',
                                              annotations={
                                                  'name': 'Generation of CO2 from nothing',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={},
                                              outflow_dict={'m01596s': [0.001, 1000]  # CO2
                                                            }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_17',
                                              annotations={
                                                  'name': 'ATP  Phosphorilated per mol of Glucose Exceeds 40 mol',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m01965s': [0, 1],  # Glucose
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={},
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1,
                                                                             'm02039c': 1},
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
                                              reaction_dict={'pim_to_pic': ({'m02751m': -1, 'm02751c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_19',
                                              annotations={
                                                  'name': 'Anaerobic Production of Propanoate from Glucose',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m01965s': [0, 1]},  # Glucose
                                              outflow_dict={'m02772s': [1, 1000]},  # Propanoate
                                              reaction_dict={}
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_20',
                                              annotations={
                                                  'name': 'Anaerobic ATP phosphorylation per glucose exceeds 2',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m01965s': [0, 1]},  # Glucose
                                              outflow_dict={},
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1,
                                                                             'm02039c': 1},
                                                                            (2.5, 1000))
                                                             }
                                              )]
    essential_tasks = essential_tasks + [Task(name='T_FAIL_21',
                                              annotations={
                                                  'name': 'ATP Production from Pi and O2',
                                                  'task_group': 'Should Fail'},
                                              should_fail=True,
                                              inflow_dict={'m02751s': [0, 1000],  # Pi
                                                           'm02630s': [0, 1000]  # O2
                                                           },
                                              outflow_dict={},
                                              reaction_dict={'atp_to_adp': ({'m01371c': -1, 'm02040c': -1,
                                                                             'm01285c': 1, 'm02751c': 1,
                                                                             'm02039c': 1},
                                                                            (1, 1000))
                                                             }
                                              )]
    jtio = JSONTaskIO()
    jtio.write_task(join(base_dir, 'Metabolic_Models/GENERAL/utility_data/metabolic_tasks_hsa_cellViability.json'),
                    essential_tasks)



    '''
    Metabolic Tasks - All tasks in Human Cells
    '''

    # other_hsa_tasks_success = []

    # Here...

    # other_hsa_tasks_fail = []

    # essential_tasks_keep = np.delete(essential_tasks, [40, 41, 42, 43, 44, 46, 47, 48, 49, 50, 51, 53, 54, 55])
    # all_hsa_tasks = np.take(essential_tasks_keep, [range(0,40), 45, 52]) + other_hsa_tasks_success\
    #                + np.take(essential_tasks_keep, [range(56,77)]) + other_hsa_tasks_fail

    # jtio = JSONTaskIO()
    # jtio.write_task(join(base_dir, 'Metabolic_Models/GENERAL/utility_data/metabolic_tasks_hsa_all.json'),
    #                all_hsa_tasks)

# from cobra.io import read_sbml_model, write_sbml_model
# from cobra.flux_analysis.variability import find_blocked_reactions
# from os.path import join
# from os import getcwd
#
## Base Directories:
# base_dir = getcwd()
#
# models_dir = join(base_dir, '0Models')
# recon3d_model_dir = join(models_dir, 'Recon3D')
#
# general_code_dir = join(base_dir, 'general_code')
# general_code_utils_dir = join(general_code_dir, 'utils')
# utils_entrez_genes_dir = join(general_code_utils_dir, 'entrez_genes')
# utils_GPRs_dir = join(general_code_utils_dir, 'GPRs')
#
# '''
# Recon3D
# '''
## Read recon3d original model:
# recon3d_original = read_sbml_model(join(recon3d_model_dir, 'Recon3D.xml.gz'))
## Get model's genes:
# gene_numbers = []
# for i in recon3d_original.genes: gene_numbers.append(i.id)
# f = open(join(utils_entrez_genes_dir, 'recon3D_genes.txt'), 'w')
# for item in gene_numbers: f.write("%s\n" % item)
# f.close()
## Get model's GPRs:
# reaction_ids = []
# for i in recon3d_original.reactions: reaction_ids.append(i.id)
# f = open(join(utils_GPRs_dir, 'recon3D_GPR.txt'), 'w')
# for item in reaction_ids: f.write("%s\t%s\n" % (item, recon3d_original.reactions.get_by_id(item).gene_reaction_rule))
# f.close()


# '''
# Recon3D_consistent
# '''
# recon3d_consistent = recon3d_original.copy()
## Remove sink (sink_) and drug module (DM_) reactions:
# recon3d_consistent.remove_reactions(
#    [r for r in recon3d_consistent.reactions if r.id[:3] == 'DM_' or r.id[:5] == 'sink_'], remove_orphans=True)
## Remove reactions that became blocked after the removal of these reactions:
# blocked_reactions = find_blocked_reactions(recon3d_consistent)
# recon3d_consistent.remove_reactions(blocked_reactions, remove_orphans=True)
## Save consistent version of Recon3D into a file:
# write_sbml_model(recon3d_consistent, join(recon3d_model_dir, 'Recon3D_consistent.xml.gz'))
## Get model's genes:
# gene_numbers = []
# for i in recon3d_consistent.genes: gene_numbers.append(i.id)
# f = open(join(utils_entrez_genes_dir, 'recon3D_consistent_genes.txt'), 'w')
# for item in gene_numbers: f.write("%s\n" % item)
# f.close()
## Get model's GPRs:
# reaction_ids = []
# for i in recon3d_consistent.reactions: reaction_ids.append(i.id)
# f = open(join(utils_GPRs_dir, 'recon3D_consistent_GPR.txt'), 'w')
# for item in reaction_ids: f.write("%s\t%s\n" % (item, recon3d_consistent.reactions.get_by_id(item).gene_reaction_rule))
# f.close()
#
#
#
# '''
# Recon3D_forTcells
# '''
# from cobra import Reaction
# recon3d_for_Tcells = recon3d_consistent.copy()
## Create a reaction for the biomass reaction from the macrophage model iAB-AM0-1410:
# new_reaction = Reaction(id='biomass_mac', name='Biomass reaction from macrophage model iAB-AM0-1410',
#                        subsystem='Exchange/demand reaction',
#                        lower_bound=0., upper_bound=1000.)
# recon3d_for_Tcells.add_reactions([new_reaction])
# recon3d_for_Tcells.reactions.biomass_mac.add_metabolites({
#    'ala_L[c]': -0.396559456, 'alpa_hs[c]': -0.011499127, 'amp[c]': -0.048664064, 'arg_L[c]': -0.325724532,
#    'asn_L[c]': -0.215407845, 'asp_L[c]': -0.282759085, 'atp[c]': -25.17352552,
#    'chsterol[c]': -0.020930954, 'cmp[c]': -0.042373167, 'cys_L[c]': -0.127154496,
#    'dag_hs[c]': -0.0036682, 'damp[c]': -0.021495345, 'dcmp[c]': -0.014937443, 'dgmp[c]': -0.014937443,
#    'dtmp[c]': -0.021495345,
#    'gln_L[c]': -0.280436629, 'glu_L[c]': -0.424428935, 'gly[c]': -0.366948135, 'glygn1[c]': -0.528027894,
#    'gmp[c]': -0.043710887,
#    'h2o[c]': -25.17352552, 'hdca[c]': -0.004850777, 'hdcea[c]': -0.001222285, 'his_L[c]': -0.153862747,
#    'ile_L[c]': -0.25953452,
#    'leu_L[c]': -0.580614138, 'lys_L[c]': -0.351852168,
#    'met_L[c]': -0.126573882,
#    'ocdca[c]': -0.004736708, 'ocdcea[c]': -0.003853116,
#    'pail_hs[c]': -0.003741686, 'pchol_hs[c]': -0.031527146, 'pe_hs[c]': -0.021107135, 'pglyc_hs[c]': -0.008918017,
#    'phe_L[c]': -0.214246617, 'pro_L[c]': -0.346626641, 'ps_hs[c]': -0.001024655,
#    'ser_L[c]': -0.476684207, 'sphmyln_hs[c]': -0.007049706,
#    'tag_hs[c]': -0.002742439, 'thr_L[c]': -0.303661194, 'trp_L[c]': -0.069673697, 'ttdca[c]': -0.00136164,
#    'tyr_L[c]': -0.156185203,
#    'ump[c]': -0.04602478,
#    'val_L[c]': -0.347207255,
#    'adp[c]': 25.17352552,
#    'h[c]': 25.17352552,
#    'pi[c]': 25.17352552
# })
# recon3d_for_Tcells.reactions.biomass_mac.gene_reaction_rule = ''
## Save this version of Recon3D into a file:
# write_sbml_model(recon3d_consistent, join(recon3d_model_dir, 'Recon3D_forTcells.xml.gz'))
## Get model's genes:
# gene_numbers = []
# for i in recon3d_consistent.genes: gene_numbers.append(i.id)
# f = open(join(utils_entrez_genes_dir, 'recon3D_forTcells_genes.txt'), 'w')
# for item in gene_numbers: f.write("%s\n" % item)
# f.close()
## Get model's GPRs:
# reaction_ids = []
# for i in recon3d_consistent.reactions: reaction_ids.append(i.id)
# f = open(join(utils_GPRs_dir, 'recon3D_forTcells_GPR.txt'), 'w')
# for item in reaction_ids: f.write("%s\t%s\n" % (item, recon3d_consistent.reactions.get_by_id(item).gene_reaction_rule))
# f.close()
