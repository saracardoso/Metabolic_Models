if __name__ == '__main__':
    from os.path import join
    import numpy as np

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
                                                  'name': 'Aerobic Rephosphorylation of GTP from Glucose',
                                                  'task_group': 'Rephosphorylation of Nucleoside Triphosphates'},
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
                                                  'name': 'Aerobic Rephosphorylation of CTP from Glucose',
                                                  'task_group': 'Rephosphorylation of Nucleoside Triphosphates'},
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
                                                  'name': 'Aerobic Rephosphorylation of UTP from Glucose',
                                                  'task_group': 'Rephosphorylation of Nucleoside Triphosphates'},
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
    essential_tasks = essential_tasks + [Task(name='T_SUCCESS_13',
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
    jtio.write_task(join(base_dir, 'GENERAL/utility_data/metabolic_tasks_hsa_cellViability.json'),
                    essential_tasks)

    '''
    Metabolic Tasks - All tasks in Human Cells
    '''

    other_hsa_tasks_success = []

    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_57',
                                                              annotations={
                                                                  'name': 'Aerobic rephosphorylation of ATP '
                                                                          'from a fatty acid',
                                                                  'task_group': 'Rephosphorilation of Nucleoside '
                                                                                'Triphosphates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02674s': [0, 1000]  # Palmitate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02039c': [0, 1000]  # H+
                                                                            },
                                                              reaction_dict={'atp_to_adp': ({'m01371c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm01285c': 1,
                                                                                             'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (1, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_58',
                                                              annotations={
                                                                  'name': 'Anaerobic rephosphorylation of ATP',
                                                                  'task_group': 'Rephosphorilation of Nucleoside '
                                                                                'Triphosphates'},
                                                              should_fail=False,
                                                              inflow_dict={'m01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm02403s': [0, 1000],  # L-Lactate
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02039c': [0, 1000]  # H+
                                                                            },
                                                              reaction_dict={'atp_to_adp': ({'m01371c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm01285c': 1, 'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (1, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_59',
                                                              annotations={
                                                                  'name': 'Anaerobic rephosphorylation of GTP',
                                                                  'task_group': 'Rephosphorilation of Nucleoside '
                                                                                'Triphosphates'},
                                                              should_fail=False,
                                                              inflow_dict={'m01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm02403s': [0, 1000],  # L-Lactate
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02039c': [0, 1000]  # H+
                                                                            },
                                                              reaction_dict={'gtp_to_gdp': ({'m02034c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm01948c': 1,
                                                                                             'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (1, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_60',
                                                              annotations={
                                                                  'name': 'Anaerobic rephosphorylation of CTP',
                                                                  'task_group': 'Rephosphorilation of Nucleoside '
                                                                                'Triphosphates'},
                                                              should_fail=False,
                                                              inflow_dict={'m01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm02403s': [0, 1000],  # L-Lactate
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02039c': [0, 1000]  # H+
                                                                            },
                                                              reaction_dict={'ctp_to_cdp': ({'m01623c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm01424c': 1, 'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (1, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_61',
                                                              annotations={
                                                                  'name': 'Anaerobic rephosphorylation of UTP',
                                                                  'task_group': 'Rephosphorilation of Nucleoside '
                                                                                'Triphosphates'},
                                                              should_fail=False,
                                                              inflow_dict={'m01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm02403s': [0, 1000],  # L-Lactate
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02039c': [0, 1000]  # H+
                                                                            },
                                                              reaction_dict={'utp_to_udp': ({'m03130c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm03106c': 1, 'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (1, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_62',
                                                              annotations={
                                                                  'name': 'ATP salvage from Adenosine',
                                                                  'task_group': 'Salvage of Nucleotides'},
                                                              should_fail=False,
                                                              inflow_dict={'m02751s': [0, 1000],  # Pi
                                                                           'm02039c': [0, 1000],  # H+
                                                                           'm01280s': [1, 1]  # Adenosine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01371s': [1, 1000]  # ATP
                                                                            },
                                                              reaction_dict={'atp_to_adp': ({'m01371c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm01285c': 1, 'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (-1000, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_63',
                                                              annotations={
                                                                  'name': 'ATP salvage from Hypoxanthine',
                                                                  'task_group': 'Salvage of Nucleotides'},
                                                              should_fail=False,
                                                              inflow_dict={'m02751s': [0, 1000],  # Pi
                                                                           'm02630s': [0, 1000],  # O2
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm02039c': [0, 1000],  # H+
                                                                           'm02159s': [1, 1],  # Hypoxanthine
                                                                           'm02806c': [1, 1]  # PRPP
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm02751s': [0, 1000],  # Pi
                                                                            'm01371s': [1, 1000]  # ATP
                                                                            },
                                                              reaction_dict={'atp_to_adp': ({'m01371c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm01285c': 1, 'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (-1000, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_64',
                                                              annotations={
                                                                  'name': 'dTTP salvage from Thymine',
                                                                  'task_group': 'Salvage of Nucleotides'},
                                                              should_fail=False,
                                                              inflow_dict={'m02751s': [0, 1000],  # Pi
                                                                           'm02630s': [0, 1000],  # O2
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm02039c': [0, 1000],  # H+
                                                                           'm02997s': [1, 1],  # Thymine
                                                                           'm00639c': [1, 1]  # 2-deoxy-D-ribose-
                                                                           # 1-phosphate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm02759s': [0, 1000],  # PPi
                                                                            'm01753c': [1, 1000]  # dTTP
                                                                            },
                                                              reaction_dict={'atp_to_adp': ({'m01371c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm01285c': 1, 'm02751c': 1,
                                                                                             'm02039c': 1},
                                                                                            (-1000, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_65',
                                                              annotations={
                                                                  'name': 'Aerobic reduction of NAD+',
                                                                  'task_group': 'Aerobic Reductions'},
                                                              should_fail=False,
                                                              inflow_dict={'m01965s': [0, 1000],  # Glucose
                                                                           'm02630s': [0, 1000]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            },
                                                              reaction_dict={'nadh_to_nad': ({'m02553c': -1,
                                                                                              'm02552c': 1},
                                                                                             (1, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_66',
                                                              annotations={
                                                                  'name': 'Aerobic reduction of NADP+',
                                                                  'task_group': 'Aerobic Reductions'},
                                                              should_fail=False,
                                                              inflow_dict={'m01965s': [0, 1000],  # Glucose
                                                                           'm02630s': [0, 1000]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            },
                                                              reaction_dict={'nadph_to_nadp': ({'m02555c': -1,
                                                                                                'm02554c': 1},
                                                                                               (1, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_67',
                                                              annotations={
                                                                  'name': 'Gluconeogenesis from Lactate',
                                                                  'task_group': 'Gluconeogenesis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm02403s': [0, 1000]  # Lactate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01965s': [1, 1000]  # Glucose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_68',
                                                              annotations={
                                                                  'name': 'Gluconeogenesis from Glycerol',
                                                                  'task_group': 'Gluconeogenesis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm01983s': [0, 1000]  # Glycerol
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01965s': [1, 1000]  # Glucose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_69',
                                                              annotations={
                                                                  'name': 'Gluconeogenesis from Alanine',
                                                                  'task_group': 'Gluconeogenesis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm01307s': [0, 1000]  # Alanine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [0, 1000],  # Urea
                                                                            'm01965s': [1, 1000]  # Glucose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_70',
                                                              annotations={
                                                                  'name': 'Gluconeogenesis from Lactate and '
                                                                          'Optionally Fatty Acid',
                                                                  'task_group': 'Gluconeogenesis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm02674s': [0, 1000],  # Palmitate
                                                                           'm02403s': [0, 1000]  # Lactate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01965s': [1, 1000]  # Glucose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_71',
                                                              annotations={
                                                                  'name': 'Gluconeogenesis from Glycerol and '
                                                                          'Optionally Fatty Acid',
                                                                  'task_group': 'Gluconeogenesis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm02674s': [0, 1000],  # Palmitate
                                                                           'm01983s': [0, 1000]  # Glycerol
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01965s': [1, 1000]  # Glucose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_72',
                                                              annotations={
                                                                  'name': 'Gluconeogenesis from Alanine and '
                                                                          'Optionally Fatty Acid',
                                                                  'task_group': 'Gluconeogenesis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02040s': [0, 1000],  # H2O
                                                                           'm02674s': [0, 1000],  # Palmitate
                                                                           'm01307s': [0, 1000]  # Alanine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01965s': [1, 1000]  # Glucose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_73',
                                                              annotations={
                                                                  'name': 'Glucose Storage in Glycogen',
                                                                  'task_group': 'Glucose Storage'},
                                                              should_fail=False,
                                                              inflow_dict={'m01965s': [11, 11],  # Glucose
                                                                           'm01996c': [1, 1],  # Glycogenin
                                                                           'm02630s': [0, 1000]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01990c': [1, 1000]  # glycogenin G11
                                                                            },
                                                              reaction_dict={'atp_to_adp': ({'m01371c': -1,
                                                                                             'm02040c': -1,
                                                                                             'm03106c': 1, 'm02751c': 1,
                                                                                             'm01285c': 1},
                                                                                            (-1000, 1000)),
                                                                             'nadh_to_nad': ({'m02553c': -1,
                                                                                              'm02552c': 1},
                                                                                             (-1000, 1000)),
                                                                             'nadph_to_nadp': ({'m02555c': -1,
                                                                                                'm02554c': 1},
                                                                                               (-1000, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_74',
                                                              annotations={
                                                                  'name': 'Release of Glucose from Glycogen',
                                                                  'task_group': 'Glucose Storage'},
                                                              should_fail=False,
                                                              inflow_dict={'m01990c': [1, 1],  # glycogenin G11
                                                                           'm02040s': [0, 1000]  # H2O
                                                                           },
                                                              outflow_dict={'m01965s': [11, 11],  # Glucose
                                                                            'm01996c': [1, 1]  # Glycogenin
                                                                            },
                                                              reaction_dict={'g1p_to_glucose': ({'m01967c': -1,
                                                                                                 'm02040c': -1,
                                                                                                 'm01965c': 1,
                                                                                                 'm02751c': 1},
                                                                                                (-1000, 1000))
                                                                             }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_75',
                                                              annotations={
                                                                  'name': 'Fructose Degradation',
                                                                  'task_group': 'Sugar Degradation'},
                                                              should_fail=False,
                                                              inflow_dict={'m01840s': [0, 1],  # Fructose
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_76',
                                                              annotations={
                                                                  'name': 'Galactose Degradation',
                                                                  'task_group': 'Sugar Degradation'},
                                                              should_fail=False,
                                                              inflow_dict={'m01910s': [0, 1],  # Galactose
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_77',
                                                              annotations={
                                                                  'name': 'UDP-galactose de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03107c': [1, 1000]  # UDP-galactose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_78',
                                                              annotations={
                                                                  'name': 'UDP-glucuronate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03109c': [1, 1000]  # UDP-glucuronate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_79',
                                                              annotations={
                                                                  'name': 'GDP-L-fucose de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01950c': [1, 1000]  # GDP-L-fucose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_80',
                                                              annotations={
                                                                  'name': 'GDP-mannose de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01951c': [1, 1000]  # GDP-mannose
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_81',
                                                              annotations={
                                                                  'name': 'UDP-N-acetyl D-galactosamine '
                                                                          'de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03110c': [1, 1000]  # UDP-N-acetyl D-
                                                                            # galactosamine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_82',
                                                              annotations={
                                                                  'name': 'CMP-N-acetylneuraminate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01592c': [1, 1000]  # CMP-N-
                                                                            # acetylneuraminate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_83',
                                                              annotations={
                                                                  'name': 'N-Acetylglucosamine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02527c': [1, 1000]  # N-Acetylglucosamine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_84',
                                                              annotations={
                                                                  'name': 'Glucuronate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01973c': [1, 1000]  # Glucuronate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_85',
                                                              annotations={
                                                                  'name': 'Alanine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01307s': [1, 1000]  # Alanine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_86',
                                                              annotations={
                                                                  'name': 'Arginine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01365s': [1, 1000]  # Arginine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_87',
                                                              annotations={
                                                                  'name': 'Asparagine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01369s': [1, 1000]  # Asparagine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_88',
                                                              annotations={
                                                                  'name': 'Aspartate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01370s': [1, 1000]  # Aspartate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_89',
                                                              annotations={
                                                                  'name': 'Glutamate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01974s': [1, 1000]  # Glutamate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_90',
                                                              annotations={
                                                                  'name': 'Glycine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01986s': [1, 1000]  # Glycine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_91',
                                                              annotations={
                                                                  'name': 'Glutamine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01975s': [1, 1000]  # Glutamine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_92',
                                                              annotations={
                                                                  'name': 'Proline de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02770s': [1, 1000]  # Proline
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_93',
                                                              annotations={
                                                                  'name': 'Serine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02896s': [1, 1000]  # Serine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_94',
                                                              annotations={
                                                                  'name': 'Cysteine de novo Synthesis from '
                                                                          'other Amino acids',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02946s': [0, 1000],  # Sulfate
                                                                           'm02125s': [0, 1000],  # Histidine
                                                                           'm02184s': [0, 1000],  # Isoleucine
                                                                           'm02360s': [0, 1000],  # Leucine
                                                                           'm02426s': [0, 1000],  # Lysine
                                                                           'm02471s': [0, 1000],  # Methionine
                                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                                           'm02993s': [0, 1000],  # Threonine
                                                                           'm03089s': [0, 1000]  # Tryptophan
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01628s': [1, 1000]  # Cysteine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_95',
                                                              annotations={
                                                                  'name': 'Tyrosine de novo Synthesis from '
                                                                          'other Amino acids',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02946s': [0, 1000],  # Sulfate
                                                                           'm02125s': [0, 1000],  # Histidine
                                                                           'm02184s': [0, 1000],  # Isoleucine
                                                                           'm02360s': [0, 1000],  # Leucine
                                                                           'm02426s': [0, 1000],  # Lysine
                                                                           'm02471s': [0, 1000],  # Methionine
                                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                                           'm02993s': [0, 1000],  # Threonine
                                                                           'm03089s': [0, 1000]  # Tryptophan
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03101s': [1, 1000]  # Tyrosine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_96',
                                                              annotations={
                                                                  'name': 'Homocysteine de novo Synthesis '
                                                                          'from other Amino acids',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02946s': [0, 1000],  # Sulfate
                                                                           'm02125s': [0, 1000],  # Histidine
                                                                           'm02184s': [0, 1000],  # Isoleucine
                                                                           'm02360s': [0, 1000],  # Leucine
                                                                           'm02426s': [0, 1000],  # Lysine
                                                                           'm02471s': [0, 1000],  # Methionine
                                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                                           'm02993s': [0, 1000],  # Threonine
                                                                           'm03089s': [0, 1000]  # Tryptophan
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02133c': [1, 1000]  # Homocystein
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_97',
                                                              annotations={
                                                                  'name': 'beta-Alanine de novo Synthesis '
                                                                          'from other Amino acids',
                                                                  'task_group': 'De novo Synthesis of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01383s': [1, 1000]  # beta-Alanine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_98',
                                                              annotations={
                                                                  'name': 'Alanine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01307s': [0, 1],  # Alanine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_99',
                                                              annotations={
                                                                  'name': 'Arginine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01365s': [0, 1],  # Arginine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_100',
                                                              annotations={
                                                                  'name': 'Asparagine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01369s': [0, 1],  # Asparagine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_101',
                                                              annotations={
                                                                  'name': 'Aspartate Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01370s': [0, 1],  # Aspartate
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_102',
                                                              annotations={
                                                                  'name': 'Cysteine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01628s': [0, 1],  # Cysteine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02042s': [0, 1000]  # H2S
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_103',
                                                              annotations={
                                                                  'name': 'Glutamate Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01974s': [0, 1],  # Glutamate
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_104',
                                                              annotations={
                                                                  'name': 'Glycine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01986s': [0, 1],  # Glycine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_105',
                                                              annotations={
                                                                  'name': 'Histidine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02125s': [0, 1],  # Histidine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_106',
                                                              annotations={
                                                                  'name': 'Isoleucine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02184s': [0, 1],  # Isoleucine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_107',
                                                              annotations={
                                                                  'name': 'Glutamine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01975s': [0, 1],  # Glutamine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_108',
                                                              annotations={
                                                                  'name': 'Leucine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02360s': [0, 1],  # Leucine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_109',
                                                              annotations={
                                                                  'name': 'Lysine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02426s': [0, 1],  # Lysine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_110',
                                                              annotations={
                                                                  'name': 'Methionine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02471s': [0, 1],  # Methionine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02042s': [0, 1000]  # H2S
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_111',
                                                              annotations={
                                                                  'name': 'Phenylalanine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02724s': [0, 1],  # Phenylalanine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_112',
                                                              annotations={
                                                                  'name': 'Proline Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02770s': [0, 1],  # Proline
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_113',
                                                              annotations={
                                                                  'name': 'Serine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02896s': [0, 1],  # Serine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_114',
                                                              annotations={
                                                                  'name': 'Threonine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02993s': [0, 1],  # Threonine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_115',
                                                              annotations={
                                                                  'name': 'Tryptophan Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m03089s': [0, 1],  # Tryptophan
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_116',
                                                              annotations={
                                                                  'name': 'Tyrosine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m03101s': [0, 1],  # Tyrosine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_117',
                                                              annotations={
                                                                  'name': 'Valine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m03135s': [0, 1],  # Valine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_118',
                                                              annotations={
                                                                  'name': 'Homocysteine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02133c': [0, 1],  # Homocysteine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02042s': [0, 1000]  # H2S
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_119',
                                                              annotations={
                                                                  'name': 'Ornithine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m02658s': [0, 1],  # Ornithine
                                                                           'm02630s': [0, 1000]  # O2
                                                                           },
                                                              outflow_dict={'m03121s': [1, 1000],  # Urea
                                                                            'm02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_120',
                                                              annotations={
                                                                  'name': 'beta-Alanine Degradation',
                                                                  'task_group': 'Degradation of Amino acis'},
                                                              should_fail=False,
                                                              inflow_dict={'m01383s': [0, 1],  # beta-Alanine
                                                                           'm02630s': [0, 1]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000]  # Urate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_121',
                                                              annotations={
                                                                  'name': 'Urea from Alanine',
                                                                  'task_group': 'Synthesis of Urea'},
                                                              should_fail=False,
                                                              inflow_dict={'m01307s': [0, 2],  # Alanine
                                                                           'm02630s': [0, 1000]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [1, 1000]  # Urea
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_122',
                                                              annotations={
                                                                  'name': 'Urea from Glutamine',
                                                                  'task_group': 'Synthesis of Urea'},
                                                              should_fail=False,
                                                              inflow_dict={'m01975s': [0, 1],  # Glutamine
                                                                           'm02630s': [0, 1000]  # O2
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [1, 1000]  # Urea
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_123',
                                                              annotations={
                                                                  'name': 'Creatine de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01619s': [1, 1000]  # Creatine
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_124',
                                                              annotations={
                                                                  'name': 'Heme de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm01821s': [0, 1000]  # Fe2+
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02049s': [1, 1000]  # Heme
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_125',
                                                              annotations={
                                                                  'name': 'PC de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm01513s': [0, 1000],  # Choline
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02684c': [1, 1000]  # PC-LD pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_126',
                                                              annotations={
                                                                  'name': 'PE de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm01797c': [0, 1000]  # Etanolamine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02685c': [1, 1000]  # PE-LD pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_127',
                                                              annotations={
                                                                  'name': 'PS de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm01513s': [0, 1000]  # Choline
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02808c': [1, 1000]  # PS-LD pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_128',
                                                              annotations={
                                                                  'name': 'PI de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm02171s': [0, 1000]  # Inositol
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02750c': [1, 1000]  # PI pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_129',
                                                              annotations={
                                                                  'name': 'Cardiolipin de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm01513s': [0, 1000]  # Choline
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01589c': [1, 1000]  # CL pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_130',
                                                              annotations={
                                                                  'name': 'SM de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Phospholipids'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm01513s': [0, 1000]  # Choline
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02908c': [1, 1000]  # SM pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_131',
                                                              annotations={
                                                                  'name': 'Ceramide de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm01513s': [0, 1000]  # Choline
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01430c': [1, 1000]  # Ceramide pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_132',
                                                              annotations={
                                                                  'name': 'Lactosylceramide de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000],  # Linolenate
                                                                           'm01513s': [0, 1000]  # Choline
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02328c': [1, 1000]  # LacCer pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_133',
                                                              annotations={
                                                                  'name': 'CoA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02471s': [0, 1000],  # Methionine
                                                                           'm02680s': [0, 1000],  # Pantothenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm01597c': [1, 1000]  # CoA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_134',
                                                              annotations={
                                                                  'name': 'NAD de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02583s': [0, 1000],  # Nicotinamide
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02553c': [1, 1000]  # NAD+
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_135',
                                                              annotations={
                                                                  'name': 'NADP de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02583s': [0, 1000]  # Nicotinamide
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02555c': [1, 1000]  # NADP+
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_136',
                                                              annotations={
                                                                  'name': 'FAD de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02842s': [0, 1000]  # Riboflavin
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01802c': [1, 1000]  # FAD
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_137',
                                                              annotations={
                                                                  'name': 'Thioredoxin de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02125s': [0, 1000],  # Histidine
                                                                           'm02184s': [0, 1000],  # Isoleucine
                                                                           'm02360s': [0, 1000],  # Leucine
                                                                           'm02426s': [0, 1000],  # Lysine
                                                                           'm02471s': [0, 1000],  # Methionine
                                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                                           'm02993s': [0, 1000],  # Threonine
                                                                           'm03089s': [0, 1000],  # Tryptophan
                                                                           'm03135s': [0, 1000]  # Valine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02042s': [0, 1000],  # H2S
                                                                            'm02990c': [1, 1000]  # Thioredoxin
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_138',
                                                              annotations={
                                                                  'name': 'THF de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm01830s': [0, 1000]  # Folate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02980c': [1, 1000]  # THF
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_139',
                                                              annotations={
                                                                  'name': 'Pyridoxal-P de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02817s': [0, 1000]  # Pyridoxine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02814c': [1, 1000]  # Pyridoxal-P 
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_140',
                                                              annotations={
                                                                  'name': 'Acetoacetate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of '
                                                                                'Other Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01253m': [1, 1000]  # Acetoacetate 
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_141',
                                                              annotations={
                                                                  'name': '(R)-3-Hydroxybutanoate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00157c': [1, 1000]
                                                                            # (R)-3-Hydroxybutanoate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_142',
                                                              annotations={
                                                                  'name': 'Farnesyl-PP de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                                           'm02751s': [0, 1000]  # Pi
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01806c': [1, 1000]  # Farnesyl-PP 
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_143',
                                                              annotations={
                                                                  'name': 'Lauric acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02344c': [1, 1000]  # Lauric acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_144',
                                                              annotations={
                                                                  'name': 'Tridecylic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03051c': [1, 1000]  # Tridecylic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_145',
                                                              annotations={
                                                                  'name': 'Myristic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02494c': [1, 1000],  # Myristic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_146',
                                                              annotations={
                                                                  'name': '9E-tetradecenoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00128c': [1, 1000],
                                                                            # 9E-tetradecenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_147',
                                                              annotations={
                                                                  'name': '7Z-tetradecenoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00117c': [1, 1000]
                                                                            # 7Z-tetradecenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_148',
                                                              annotations={
                                                                  'name': 'Physeteric acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02745c': [1, 1000]  # Physeteric acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_149',
                                                              annotations={
                                                                  'name': 'Pentadecylic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02690c': [1, 1000]  # Pentadecylic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_150',
                                                              annotations={
                                                                  'name': 'Palmitate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02674c': [1, 1000]  # Palmitate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_151',
                                                              annotations={
                                                                  'name': 'Palmitolate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02675c': [1, 1000]  # Palmitolate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_152',
                                                              annotations={
                                                                  'name': '7-palmitoleic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01197c': [1, 1000]  # 7-palmitoleic acid 
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_153',
                                                              annotations={
                                                                  'name': 'Margaric acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02456c': [1, 1000]  # Margaric acid 
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_154',
                                                              annotations={
                                                                  'name': '10Z-heptadecenoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm00003c': [1, 1000]
                                                                            # 10Z-heptadecenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_155',
                                                              annotations={
                                                                  'name': '9-heptadecylenic  acid  de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm01238c': [1, 1000]
                                                                            # 9-heptadecylenic  acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_156',
                                                              annotations={
                                                                  'name': 'Stearate  de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02938c': [1, 1000]  # Stearate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_157',
                                                              annotations={
                                                                  'name': '13Z-octadecenoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00019c': [1, 1000]
                                                                            # 13Z-octadecenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_158',
                                                              annotations={
                                                                  'name': 'cis-vaccenic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01585c': [1, 1000]  # cis-vaccenic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_159',
                                                              annotations={
                                                                  'name': 'Oleate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02646c': [1, 1000]  # Oleate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_160',
                                                              annotations={
                                                                  'name': 'Elaidate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01778c': [1, 1000]  # Elaidate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_161',
                                                              annotations={
                                                                  'name': '7Z-octadecenoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00115c': [1, 1000]  # 7Z-octadecenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_162',
                                                              annotations={
                                                                  'name': '6Z,9Z-octadecadienoic acid de novo '
                                                                          'Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00104c': [1, 1000]
                                                                            # 6Z,9Z-octadecadienoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_163',
                                                              annotations={
                                                                  'name': 'Nonadecylic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00017c': [1, 1000]  # Nonadecylic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_164',
                                                              annotations={
                                                                  'name': 'Eicosanoate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01771c': [1, 1000]  # Eicosanoate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_165',
                                                              annotations={
                                                                  'name': '13Z-Eicosenoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00017c': [1, 1000]  # 13Z-Eicosenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_166',
                                                              annotations={
                                                                  'name': 'cis-gondoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01584c': [1, 1000]  # cis-gondoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_167',
                                                              annotations={
                                                                  'name': '9-Eicosenoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01235c': [1, 1000]  # 9-Eicosenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_168',
                                                              annotations={
                                                                  'name': '8,11-Eicosadienoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01207c': [1, 1000]
                                                                            # 8,11-Eicosadienoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_169',
                                                              annotations={
                                                                  'name': 'Mead acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02457c': [1, 1000]  # Mead acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_170',
                                                              annotations={
                                                                  'name': 'Henicosanoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02053c': [1, 1000]  # Henicosanoic acid 
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_171',
                                                              annotations={
                                                                  'name': 'Behenic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01373c': [1, 1000]  # Behenic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_172',
                                                              annotations={
                                                                  'name': 'cis-erucic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01583c': [1, 1000]  # cis-erucic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_173',
                                                              annotations={
                                                                  'name': 'cis-cetoleic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01582c': [1, 1000]  # cis-cetoleic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_174',
                                                              annotations={
                                                                  'name': 'Tricosanoic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02184c': [0, 1000]  # Isoleucine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm03045c': [1, 1000]  # Tricosanoic acid 
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_175',
                                                              annotations={
                                                                  'name': 'Lignocerate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02385c': [1, 1000]  # Lignocerate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_176',
                                                              annotations={
                                                                  'name': 'Nervonic de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02564c': [1, 1000]  # Nervonic
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_177',
                                                              annotations={
                                                                  'name': 'Cerotic de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01432c': [1, 1000]  # Cerotic
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_178',
                                                              annotations={
                                                                  'name': 'Ximenic de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000]  # Glucose
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03153c': [1, 1000]  # Ximenic
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_179',
                                                              annotations={
                                                                  'name': 'Stearidonic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02939c': [1, 1000]  # Stearidonic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_180',
                                                              annotations={
                                                                  'name': 'omega-3-Arachidonic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02648c': [1, 1000]
                                                                            # omega-3-Arachidonic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_181',
                                                              annotations={
                                                                  'name': 'EPA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01784c': [1, 1000]  # EPA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_182',
                                                              annotations={
                                                                  'name': 'DPA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01741c': [1, 1000]  # DPA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_183',
                                                              annotations={
                                                                  'name': '9Z,12Z,15Z,18Z,21Z-TPA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00135c': [1, 1000]
                                                                            # 9Z,12Z,15Z,18Z,21Z-TPA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_184',
                                                              annotations={
                                                                  'name': '6Z,9Z,12Z,15Z,18Z,21Z-THA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00114c': [1, 1000]
                                                                            # 6Z,9Z,12Z,15Z,18Z,21Z-THA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_185',
                                                              annotations={
                                                                  'name': 'DHA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01689c': [1, 1000]  # DHA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_186',
                                                              annotations={
                                                                  'name': '11Z,14Z,17Z-eicosatrienoic acid de '
                                                                          'novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00010c': [1, 1000]
                                                                            # 11Z,14Z,17Z-eicosatrienoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_187',
                                                              annotations={
                                                                  'name': '13,16,19-docosatrienoic acid de '
                                                                          'novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00341c': [1, 1000]
                                                                            # 13,16,19-docosatrienoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_188',
                                                              annotations={
                                                                  'name': '10,13,16,19-docosatetraenoic acid de '
                                                                          'novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00260c': [1, 1000]
                                                                            # 10,13,16,19-docosatetraenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_189',
                                                              annotations={
                                                                  'name': '12,15,18,21-tetracosatetraenoic acid de '
                                                                          'novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00315c': [1, 1000]
                                                                            # 12,15,18,21-tetracosatetraenoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_190',
                                                              annotations={
                                                                  'name': 'gamma-Linolenate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01932c': [1, 1000]  # gamma-Linolenate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_191',
                                                              annotations={
                                                                  'name': 'Dihomo-gamma-linolenate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01696c': [1, 1000]
                                                                            # Dihomo-gamma-linolenate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_192',
                                                              annotations={
                                                                  'name': 'Arachidonate de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01362c': [1, 1000]  # Arachidonate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_193',
                                                              annotations={
                                                                  'name': 'Adrenic acid de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01291c': [1, 1000]  # Adrenic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_194',
                                                              annotations={
                                                                  'name': '9Z,12Z,15Z,18Z-TTA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00132c': [1, 1000]  # 9Z,12Z,15Z,18Z-TTA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_195',
                                                              annotations={
                                                                  'name': '6Z,9Z,12Z,15Z,18Z-TPA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00111c': [1, 1000]
                                                                            # 6Z,9Z,12Z,15Z,18Z-TPA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_196',
                                                              annotations={
                                                                  'name': '4Z,7Z,10Z,13Z,16Z-DPA de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00094c': [1, 1000]
                                                                            # 4Z,7Z,10Z,13Z,16Z-DPA
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_197',
                                                              annotations={
                                                                  'name': '11Z,14Z-eicosadienoic acid de '
                                                                          'novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00008c': [1, 1000]
                                                                            # 11Z,14Z-eicosadienoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_198',
                                                              annotations={
                                                                  'name': '13Z,16Z-docosadienoic acid de '
                                                                          'novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00021c': [1, 1000]
                                                                            # 13Z,16Z-docosadienoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_199',
                                                              annotations={
                                                                  'name': '10,13,16-docosatriynoic acid de '
                                                                          'novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Other '
                                                                                'Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm00265c': [1, 1000]
                                                                            # 10,13,16-docosatriynoic acid
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_200',
                                                              annotations={
                                                                  'name': 'Lauric acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02344s': [1, 1]  # Lauric acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_201',
                                                              annotations={
                                                                  'name': 'Tridecylic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm03051s': [1, 1]  # Tridecylic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_202',
                                                              annotations={
                                                                  'name': 'Myristic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02494s': [1, 1]  # Myristic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_203',
                                                              annotations={
                                                                  'name': '9E-tetradecenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00128s': [1, 1]  # 9E-tetradecenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_204',
                                                              annotations={
                                                                  'name': '7Z-tetradecenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00117s': [1, 1]  # 7Z-tetradecenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_205',
                                                              annotations={
                                                                  'name': 'Physeteric acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02745s': [1, 1]  # Physeteric acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_206',
                                                              annotations={
                                                                  'name': 'Pentadecylic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02690s': [1, 1]  # Pentadecylic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_207',
                                                              annotations={
                                                                  'name': 'Palmitate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02674s': [1, 1]  # Palmitate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_208',
                                                              annotations={
                                                                  'name': 'Palmitolate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02675s': [1, 1]  # Palmitolate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_209',
                                                              annotations={
                                                                  'name': '7-palmitoleic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01197s': [1, 1]  # 7-palmitoleic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_210',
                                                              annotations={
                                                                  'name': 'Margaric acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02456s': [1, 1]  # Margaric acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_211',
                                                              annotations={
                                                                  'name': '10Z-heptadecenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00003s': [1, 1]  # 10Z-heptadecenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_212',
                                                              annotations={
                                                                  'name': '9-heptadecylenic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01238s': [1, 1]  # 9-heptadecylenic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_213',
                                                              annotations={
                                                                  'name': 'Stearate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02938s': [1, 1]  # Stearate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_214',
                                                              annotations={
                                                                  'name': '13Z-octadecenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00019s': [1, 1]  # 13Z-octadecenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_215',
                                                              annotations={
                                                                  'name': 'cis-vaccenic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01585s': [1, 1]  # cis-vaccenic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_216',
                                                              annotations={
                                                                  'name': 'Oleate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02646s': [1, 1]  # Oleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_217',
                                                              annotations={
                                                                  'name': 'Elaidate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01778s': [1, 1]  # Elaidate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_218',
                                                              annotations={
                                                                  'name': '7Z-octadecenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00115s': [1, 1]  # 7Z-octadecenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_219',
                                                              annotations={
                                                                  'name': '6Z,9Z-octadecadienoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00104s': [1, 1]
                                                                           # 6Z,9Z-octadecadienoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_220',
                                                              annotations={
                                                                  'name': 'Nonadecylic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02613s': [1, 1]  # Nonadecylic
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_221',
                                                              annotations={
                                                                  'name': 'Nonadecylic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02613s': [1, 1]  # Nonadecylic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_222',
                                                              annotations={
                                                                  'name': 'Eicosanoate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01771s': [1, 1]  # Eicosanoate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_223',
                                                              annotations={
                                                                  'name': '13Z-Eicosenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00017s': [1, 1]  # 13Z-Eicosenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_224',
                                                              annotations={
                                                                  'name': 'cis-gondoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01584s': [1, 1]  # cis-gondoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_225',
                                                              annotations={
                                                                  'name': '9-Eicosenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01235s': [1, 1]  # 9-Eicosenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_226',
                                                              annotations={
                                                                  'name': '8,11-Eicosadienoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01207s': [1, 1]  # 8,11-Eicosadienoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_227',
                                                              annotations={
                                                                  'name': 'Mead acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02457s': [1, 1]  # Mead acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_228',
                                                              annotations={
                                                                  'name': 'Henicosanoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02053s': [1, 1]  # Henicosanoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_229',
                                                              annotations={
                                                                  'name': 'Behenic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01373s': [1, 1]  # Behenic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_230',
                                                              annotations={
                                                                  'name': 'cis-erucic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01583s': [1, 1]  # cis-erucic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_231',
                                                              annotations={
                                                                  'name': 'cis-cetoleic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01582s': [1, 1]  # cis-cetoleic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_232',
                                                              annotations={
                                                                  'name': 'Tricosanoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm03045s': [1, 1]  # Tricosanoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_233',
                                                              annotations={
                                                                  'name': 'Lignocerate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02385s': [1, 1]  # Lignocerate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_234',
                                                              annotations={
                                                                  'name': 'Nervonic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02564s': [1, 1]  # Nervonic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_235',
                                                              annotations={
                                                                  'name': 'Cerotic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01432s': [1, 1]  # Cerotic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_236',
                                                              annotations={
                                                                  'name': 'Ximenic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm03153s': [1, 1]  # Ximenic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_237',
                                                              annotations={
                                                                  'name': 'Linolenate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02389s': [1, 1]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_238',
                                                              annotations={
                                                                  'name': 'Stearidonic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02939s': [1, 1]  # Stearidonic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_239',
                                                              annotations={
                                                                  'name': 'omega-3-Arachidonic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02648s': [1, 1]  # omega-3-Arachidonic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_240',
                                                              annotations={
                                                                  'name': 'EPA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01784s': [1, 1]  # EPA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_241',
                                                              annotations={
                                                                  'name': 'DPA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01741s': [1, 1]  # DPA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_242',
                                                              annotations={
                                                                  'name': '9Z,12Z,15Z,18Z,21Z-TPA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00135s': [1, 1]  # 9Z,12Z,15Z,18Z,21Z-TPA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_243',
                                                              annotations={
                                                                  'name': '6Z,9Z,12Z,15Z,18Z,21Z-THA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00114s': [1, 1]
                                                                           # 6Z,9Z,12Z,15Z,18Z,21Z-THA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_244',
                                                              annotations={
                                                                  'name': 'DHA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01689s': [1, 1]  # DHA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_245',
                                                              annotations={
                                                                  'name': '11Z,14Z,17Z-eicosatrienoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00010s': [1, 1]
                                                                           # 11Z,14Z,17Z-eicosatrienoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_246',
                                                              annotations={
                                                                  'name': '13,16,19-docosatrienoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00341s': [1, 1]
                                                                           # 13,16,19-docosatrienoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_247',
                                                              annotations={
                                                                  'name': '10,13,16,19-docosatetraenoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00260s': [1, 1]
                                                                           # 10,13,16,19-docosatetraenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_248',
                                                              annotations={
                                                                  'name': '12,15,18,21-tetracosatetraenoic '
                                                                          'acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00315s': [1, 1]
                                                                           # 12,15,18,21-tetracosatetraenoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_249',
                                                              annotations={
                                                                  'name': 'Linoleate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02387s': [1, 1]  # Linoleate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_250',
                                                              annotations={
                                                                  'name': 'gamma-Linolenate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01932s': [1, 1]  # gamma-Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_251',
                                                              annotations={
                                                                  'name': 'Dihomo-gamma-linolenate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01696s': [1, 1]  # Dihomo-gamma-linolenate 
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_252',
                                                              annotations={
                                                                  'name': 'Arachidonate Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01362s': [1, 1]  # Arachidonate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_253',
                                                              annotations={
                                                                  'name': 'Adrenic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01291s': [1, 1]  # Adrenic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_254',
                                                              annotations={
                                                                  'name': '9Z,12Z,15Z,18Z-TTA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00132s': [1, 1]  # 9Z,12Z,15Z,18Z-TTA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_255',
                                                              annotations={
                                                                  'name': '6Z,9Z,12Z,15Z,18Z-TPA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00111s': [1, 1]  # 6Z,9Z,12Z,15Z,18Z-TPA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_256',
                                                              annotations={
                                                                  'name': '4Z,7Z,10Z,13Z,16Z-DPA Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00094s': [1, 1]  # 4Z,7Z,10Z,13Z,16Z-DPA
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_257',
                                                              annotations={
                                                                  'name': '11Z,14Z-eicosadienoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00008s': [1, 1]
                                                                           # 11Z,14Z-eicosadienoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_258',
                                                              annotations={
                                                                  'name': '13Z,16Z-docosadienoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00021s': [1, 1]
                                                                           # 13Z,16Z-docosadienoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_259',
                                                              annotations={
                                                                  'name': '10,13,16-docosatriynoic acid Oxidation',
                                                                  'task_group': 'Complete Oxidations'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm00265s': [1, 1]
                                                                           # 10,13,16-docosatriynoic acid
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_260',
                                                              annotations={
                                                                  'name': 'Triacylglycerol de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02387s': [0, 1000],  # Linoleate
                                                                           'm02389s': [0, 1000]  # Linolenate
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm02958c': [1, 1000]  # TAG-LD pool
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_261',
                                                              annotations={
                                                                  'name': 'Glycocholate de novo Synthesis and '
                                                                          'Excretion',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02724s': [0, 1000]  # Phenylalanine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [0, 1000],  # Urea
                                                                            'm01988s': [1, 1000]  # Glycocholate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_262',
                                                              annotations={
                                                                  'name': 'Glycochenodeoxycholate de novo Synthesis '
                                                                          'and Excretion',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02724s': [0, 1000]  # Phenylalanine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [0, 1000],  # Urea
                                                                            'm01987s': [1, 1000]
                                                                            # Glycochenodeoxycholate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_263',
                                                              annotations={
                                                                  'name': 'Taurocholate de novo Synthesis and '
                                                                          'Excretion',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                                           'm01628s': [0, 1000]  # Cysteine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [0, 1000],  # Urea
                                                                            'm02963s': [1, 1000]  # Taurocholate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_264',
                                                              annotations={
                                                                  'name': 'Taurochenodeoxycholate de novo Synthesis '
                                                                          'and Excretion',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02724s': [0, 1000],  # Phenylalanine
                                                                           'm01628s': [0, 1000]  # Cysteine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [0, 1000],  # Urea
                                                                            'm02962s': [1, 1000]
                                                                            # Taurochenodeoxycholate
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_265',
                                                              annotations={
                                                                  'name': 'PAPS de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm02946s': [0, 1000],  # Sulfate
                                                                           'm02125s': [0, 1000],  # Histidine
                                                                           'm03135s': [0, 1000]  # Valine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02682c': [1, 1000]  # PAPS
                                                                            }
                                                              )]

    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_266',
                                                              annotations={
                                                                  'name': 'PAP Degradation',
                                                                  'task_group': 'Degradation of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02682c': [0, 1000]  # PAPS
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02751s': [0, 1000],  # Pi
                                                                            'm02042s': [0, 1000]  # H2S
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_267',
                                                              annotations={
                                                                  'name': 'SAM de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02471s': [0, 1000]  # Methionine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02877c': [1, 1000]  # SAM
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_268',
                                                              annotations={
                                                                  'name': 'GSH de novo Synthesis',
                                                                  'task_group': 'De novo Synthesis of Intermediates'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm02578s': [0, 1000],  # NH3
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02471s': [0, 1000]  # Methionine
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03120s': [0, 1000],  # Urate
                                                                            'm02026c': [1, 1000]  # GSH
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_269',
                                                              annotations={
                                                                  'name': 'Bilirubin Conjugation',
                                                                  'task_group': 'Others'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02751s': [0, 1000],  # Pi
                                                                           'm01396s': [0, 1000]  # Bilirubin
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm01397s': [1, 1000]
                                                                            # bilirubin-bisglucuronoside
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_270',
                                                              annotations={
                                                                  'name': 'NH3 Import and Degradation',
                                                                  'task_group': 'Others'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01965s': [0, 1000],  # Glucose
                                                                           'm02578s': [0, 1000]  # HN3
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000],  # CO2
                                                                            'm03121s': [1, 1000]  # Urea
                                                                            }
                                                              )]
    other_hsa_tasks_success = other_hsa_tasks_success + [Task(name='T_SUCCESS_271',
                                                              annotations={
                                                                  'name': 'Ethanol Import and Degradation',
                                                                  'task_group': 'Others'},
                                                              should_fail=False,
                                                              inflow_dict={'m02630s': [0, 1000],  # O2
                                                                           'm01796s': [1, 1]  # Ethanol
                                                                           },
                                                              outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                            'm01596s': [0, 1000]  # CO2
                                                                            }
                                                              )]

    other_hsa_tasks_fail = []
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_22',
                                                        annotations={
                                                            'name': 'Histidine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02125c': [1, 1000]  # Histidine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_23',
                                                        annotations={
                                                            'name': 'Isoleucine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02184c': [1, 1000]  # Isoleucine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_24',
                                                        annotations={
                                                            'name': 'Leucine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02360c': [1, 1000]  # Leucine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_25',
                                                        annotations={
                                                            'name': 'Lysine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02426c': [1, 1000]  # Lysine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_26',
                                                        annotations={
                                                            'name': 'Methionine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02471c': [1, 1000]  # Methionine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_27',
                                                        annotations={
                                                            'name': 'Phenylalanine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02724c': [1, 1000]  # Phenylalanine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_28',
                                                        annotations={
                                                            'name': 'Threonine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02993c': [1, 1000]  # Threonine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_29',
                                                        annotations={
                                                            'name': 'Tryptophan de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm03089c': [1, 1000]  # Tryptophan
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_30',
                                                        annotations={
                                                            'name': 'Valine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm03135c': [1, 1000]  # Valine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_31',
                                                        annotations={
                                                            'name': 'Cysteine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000],  # Pi
                                                                     'm02946s': [0, 1000]  # Sulfate
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm01628s': [1, 1000]  # Cysteine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_32',
                                                        annotations={
                                                            'name': 'Cystine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000],  # Pi
                                                                     'm02949s': [0, 1000]  # Sulfite
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm01629s': [1, 1000]  # Cystine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_33',
                                                        annotations={
                                                            'name': 'Tyrosine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm03101s': [1, 1000]  # Cystine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_34',
                                                        annotations={
                                                            'name': 'Homocysteine de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Amino acis'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000],  # Glucose
                                                                     'm02578s': [0, 1000],  # NH3
                                                                     'm02751s': [0, 1000]  # Pi
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02133c': [1, 1000]  # Homocysteine
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_35',
                                                        annotations={
                                                            'name': 'Linolenate de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Other Intermediates'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000]  # Glucose
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02389c': [1, 1000]  # Linolenate
                                                                      }
                                                        )]
    other_hsa_tasks_fail = other_hsa_tasks_fail + [Task(name='T_FAIL_36',
                                                        annotations={
                                                            'name': 'Linoleate de novo Synthesis',
                                                            'task_group': 'De novo Synthesis of Other Intermediates'},
                                                        should_fail=True,
                                                        inflow_dict={'m02630s': [0, 1000],  # O2
                                                                     'm01965s': [0, 1000]  # Glucose
                                                                     },
                                                        outflow_dict={'m02040s': [0, 1000],  # H2O
                                                                      'm01596s': [0, 1000],  # CO2
                                                                      'm02387c': [1, 1000]  # Linoleate
                                                                      }
                                                        )]

    all_hsa_tasks = np.take(essential_tasks, [*range(0, 40), 45, 52]).tolist() + other_hsa_tasks_success \
        + np.take(essential_tasks, [*range(56, len(essential_tasks))]).tolist() + other_hsa_tasks_fail

    jtio = JSONTaskIO()
    jtio.write_task(join(base_dir, 'GENERAL/utility_data/metabolic_tasks_hsa_all.json'),
                    all_hsa_tasks)


    '''
    HumanGEM GPRs
    '''

    # Directories:
    models_dir = join(base_dir, '0MODELS')
    HumanGEM_dir = join(models_dir, 'HumanGEM')

    general_utility_data_dir = join(base_dir, 'GENERAL/utility_data')

    # --- Read HumanGEM_forTcells original model ---
    print('Reading HumanGEM-1.4.1...')
    HumanGEM = read_sbml_model(join(HumanGEM_dir, 'HumanGEM-1.4.1.xml.gz'))

    # --- Get HumanGEM model's GPRs ---
    print('Collecting GPR rules from HumanGEM-1.4.1')
    reaction_ids = []
    for i in HumanGEM.reactions:
        reaction_ids.append(i.id)
    f = open(join(general_utility_data_dir, 'HumanGEM-1.4.1_GPRs.txt'), 'w')
    for item in reaction_ids:
        f.write("%s\t%s\n" % (item, HumanGEM.reactions.get_by_id(item).gene_reaction_rule))
    f.close()
