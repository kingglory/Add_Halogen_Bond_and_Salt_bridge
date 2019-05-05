-------------------------------------------------------------------------------
Initialization, inputs
**********************
  silent = False
  write_geo_file = True
  file_name = None
  show_states = False
  restraints = None
  restraints_directory = None
  output_file_name_prefix = None
  directory = None
  job_title = None
  fix_rotamer_outliers = True
  allow_allowed_rotamers = True
  stop_for_unknowns = True
  pdb_interpretation {
    restraints_library {
      cdl = True
      omega_cdl = False
      rdl = False
      hpdl = False
    }
    sort_atoms = True
    superpose_ideal_ligand = *None all SF4 F3S DVT
    flip_symmetric_amino_acids = False
    disable_uc_volume_vs_n_atoms_check = False
    correct_hydrogens = True
    secondary_structure {
      protein {
        enabled = True
        search_method = *ksdssp mmtbx_dssp from_ca cablam
        distance_ideal_n_o = 2.9
        distance_cut_n_o = 3.5
        remove_outliers = True
        restrain_hbond_angles = True
        helix {
          serial_number = None
          helix_identifier = None
          enabled = True
          selection = None
          helix_type = *alpha pi 3_10 unknown
          sigma = 0.05
          slack = 0
          top_out = False
          hbond {
            donor = None
            acceptor = None
          }
        }
        sheet {
          enabled = True
          first_strand = None
          sheet_id = None
          strand {
            selection = None
            sense = parallel antiparallel *unknown
            bond_start_current = None
            bond_start_previous = None
          }
          sigma = 0.05
          slack = 0
          top_out = False
          hbond {
            donor = None
            acceptor = None
          }
        }
      }
      nucleic_acid {
        enabled = True
        hbond_distance_cutoff = 3.4
        angle_between_bond_and_nucleobase_cutoff = 35.0
        base_pair {
          enabled = True
          base1 = None
          base2 = None
          saenger_class = 0
          restrain_planarity = False
          planarity_sigma = 0.176
          restrain_hbonds = True
          restrain_hb_angles = True
          restrain_parallelity = True
          parallelity_target = 0
          parallelity_sigma = 0.0335
        }
        stacking_pair {
          enabled = True
          base1 = None
          base2 = None
          angle = 0
          sigma = 0.027
        }
      }
      ss_by_chain = True
      max_rmsd = 1
      use_representative_chains = True
      max_representative_chains = 100
      enabled = False
    }
    c_beta_restraints = True
    reference_coordinate_restraints {
      enabled = False
      exclude_outliers = True
      selection = all
      sigma = 0.2
      limit = 1.0
      top_out = False
    }
    automatic_linking {
      link_all = False
      link_none = False
      link_metals = False
      link_residues = False
      link_amino_acid_rna_dna = False
      link_carbohydrates = True
      link_ligands = True
      link_small_molecules = False
      metal_coordination_cutoff = 3.5
      amino_acid_bond_cutoff = 1.9
      inter_residue_bond_cutoff = 2.2
      buffer_for_second_row_elements = 0.5
      carbohydrate_bond_cutoff = 1.99
      ligand_bond_cutoff = 1.99
      small_molecule_bond_cutoff = 1.98
    }
    include_in_automatic_linking {
      selection_1 = None
      selection_2 = None
      bond_cutoff = 4.5
    }
    exclude_from_automatic_linking {
      selection_1 = None
      selection_2 = None
    }
    use_neutron_distances = False
    apply_cis_trans_specification {
      cis_trans_mod = cis *trans
      residue_selection = None
    }
    apply_cif_restraints {
      restraints_file_name = None
      residue_selection = None
    }
    apply_cif_modification {
      data_mod = None
      residue_selection = None
    }
    apply_cif_link {
      data_link = None
      residue_selection_1 = None
      residue_selection_2 = None
    }
    disulfide_bond_exclusions_selection_string = None
    exclusion_distance_cutoff = 3
    link_distance_cutoff = 3
    disulfide_distance_cutoff = 3
    add_angle_and_dihedral_restraints_for_disulfides = True
    dihedral_function_type = *determined_by_sign_of_periodicity \
                             all_sinusoidal all_harmonic
    chir_volume_esd = 0.2
    peptide_link {
      ramachandran_restraints = False
      cis_threshold = 45
      apply_all_trans = False
      discard_omega = False
      discard_psi_phi = True
      apply_peptide_plane = False
      omega_esd_override_value = None
      rama_weight = 1.0
      scale_allowed = 1.0
      rama_potential = *oldfield emsley
      oldfield {
        esd = 10.0
        weight_scale = 1.0
        dist_weight_max = 10.0
        weight = None
        plot_cutoff = 0.027
      }
      rama_selection = None
      restrain_rama_outliers = True
      restrain_rama_allowed = True
      restrain_allowed_outliers_with_emsley = False
    }
    max_reasonable_bond_distance = 50.0
    nonbonded_distance_cutoff = None
    default_vdw_distance = 1
    min_vdw_distance = 1
    nonbonded_buffer = 1
    nonbonded_weight = None
    const_shrink_donor_acceptor = 0.6
    vdw_1_4_factor = 0.8
    min_distance_sym_equiv = 0.5
    custom_nonbonded_symmetry_exclusions = None
    translate_cns_dna_rna_residue_names = None
    proceed_with_excessive_length_bonds = False
    rna_sugar_pucker_analysis {
      bond_min_distance = 1.2
      bond_max_distance = 1.8
      epsilon_range_min = 155.0
      epsilon_range_max = 310.0
      delta_range_2p_min = 129.0
      delta_range_2p_max = 162.0
      delta_range_3p_min = 65.0
      delta_range_3p_max = 104.0
      p_distance_c1p_outbound_line_2p_max = 2.9
      o3p_distance_c1p_outbound_line_2p_max = 2.4
      bond_detection_distance_tolerance = 0.5
    }
    show_histogram_slots {
      bond_lengths = 5
      nonbonded_interaction_distances = 5
      bond_angle_deviations_from_ideal = 5
      dihedral_angle_deviations_from_ideal = 5
      chiral_volume_deviations_from_ideal = 5
    }
    show_max_items {
      not_linked = 5
      bond_restraints_sorted_by_residual = 5
      nonbonded_interactions_sorted_by_model_distance = 5
      bond_angle_restraints_sorted_by_residual = 5
      dihedral_angle_restraints_sorted_by_residual = 3
      chirality_restraints_sorted_by_residual = 3
      planarity_restraints_sorted_by_residual = 3
      residues_with_excluded_nonbonded_symmetry_interactions = 12
      fatal_problem_max_lines = 10
    }
    ncs_group {
      reference = None
      selection = None
    }
    ncs_search {
      enabled = False
      exclude_selection = "element H or element D or water"
      chain_similarity_threshold = 0.85
      chain_max_rmsd = 2.
      residue_match_radius = 4.0
      try_shortcuts = False
      minimum_number_of_atoms_in_copy = 3
    }
    clash_guard {
      nonbonded_distance_threshold = 0.5
      max_number_of_distances_below_threshold = 100
      max_fraction_of_distances_below_threshold = 0.1
    }
  }
  geometry_restraints {
    edits {
      excessive_bond_distance_limit = 10
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 1 and name O
        atom_selection_2 = chain A and resseq 5 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 2 and name O
        atom_selection_2 = chain A and resseq 6 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 3 and name O
        atom_selection_2 = chain A and resseq 7 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 4 and name O
        atom_selection_2 = chain A and resseq 8 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 5 and name O
        atom_selection_2 = chain A and resseq 9 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 6 and name O
        atom_selection_2 = chain A and resseq 10 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 7 and name O
        atom_selection_2 = chain A and resseq 11 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 8 and name O
        atom_selection_2 = chain A and resseq 12 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 9 and name O
        atom_selection_2 = chain A and resseq 13 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 10 and name O
        atom_selection_2 = chain A and resseq 14 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 11 and name O
        atom_selection_2 = chain A and resseq 15 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 12 and name O
        atom_selection_2 = chain A and resseq 16 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 13 and name O
        atom_selection_2 = chain A and resseq 17 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 14 and name O
        atom_selection_2 = chain A and resseq 18 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 15 and name O
        atom_selection_2 = chain A and resseq 19 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      bond {
        action = *add delete change
        atom_selection_1 = chain A and resseq 16 and name O
        atom_selection_2 = chain A and resseq 20 and name N
        symmetry_operation = None
        distance_ideal = 2.900000
        sigma = 0.1
        slack = None
        limit = -0.1
        top_out = False
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 1 and name O
        atom_selection_2 = chain A and resseq 5 and name N
        atom_selection_3 = chain A and resseq 1 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 2 and name O
        atom_selection_2 = chain A and resseq 6 and name N
        atom_selection_3 = chain A and resseq 2 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 3 and name O
        atom_selection_2 = chain A and resseq 7 and name N
        atom_selection_3 = chain A and resseq 3 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 4 and name O
        atom_selection_2 = chain A and resseq 8 and name N
        atom_selection_3 = chain A and resseq 4 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 5 and name O
        atom_selection_2 = chain A and resseq 9 and name N
        atom_selection_3 = chain A and resseq 5 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 6 and name O
        atom_selection_2 = chain A and resseq 10 and name N
        atom_selection_3 = chain A and resseq 6 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 7 and name O
        atom_selection_2 = chain A and resseq 11 and name N
        atom_selection_3 = chain A and resseq 7 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 8 and name O
        atom_selection_2 = chain A and resseq 12 and name N
        atom_selection_3 = chain A and resseq 8 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 9 and name O
        atom_selection_2 = chain A and resseq 13 and name N
        atom_selection_3 = chain A and resseq 9 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 10 and name O
        atom_selection_2 = chain A and resseq 14 and name N
        atom_selection_3 = chain A and resseq 10 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 11 and name O
        atom_selection_2 = chain A and resseq 15 and name N
        atom_selection_3 = chain A and resseq 11 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 12 and name O
        atom_selection_2 = chain A and resseq 16 and name N
        atom_selection_3 = chain A and resseq 12 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 13 and name O
        atom_selection_2 = chain A and resseq 17 and name N
        atom_selection_3 = chain A and resseq 13 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 14 and name O
        atom_selection_2 = chain A and resseq 18 and name N
        atom_selection_3 = chain A and resseq 14 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 15 and name O
        atom_selection_2 = chain A and resseq 19 and name N
        atom_selection_3 = chain A and resseq 15 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      angle {
        action = *add delete change
        atom_selection_1 = chain A and resseq 16 and name O
        atom_selection_2 = chain A and resseq 20 and name N
        atom_selection_3 = chain A and resseq 16 and name C
        angle_ideal = 147.150000
        sigma = 5
      }
      dihedral {
        action = *add delete change
        atom_selection_1 = None
        atom_selection_2 = None
        atom_selection_3 = None
        atom_selection_4 = None
        angle_ideal = None
        sigma = None
        periodicity = None
      }
      planarity {
        action = *add delete change
        atom_selection = None
        sigma = None
      }
      parallelity {
        action = *add delete change
        atom_selection_1 = None
        atom_selection_2 = None
        sigma = 0.027
        target_angle_deg = 0
      }
    }
    remove {
      angles = None
      dihedrals = None
      chiralities = None
      planarities = None
      parallelities = None
    }
  }
  reference_model {
    enabled = False
    file = None
    use_starting_model_as_reference = False
    sigma = 1.0
    limit = 15.0
    hydrogens = False
    main_chain = True
    side_chain = True
    fix_outliers = True
    strict_rotamer_matching = False
    auto_shutoff_for_ncs = False
    secondary_structure_only = False
    reference_group {
      reference = None
      selection = None
      file_name = None
    }
    search_options {
      exclude_selection = "element H or element D or water"
      chain_similarity_threshold = 0.85
      chain_max_rmsd = 100.
      residue_match_radius = 1000
      try_shortcuts = False
      minimum_number_of_atoms_in_copy = 3
    }
  }
  selection = all
  minimization {
    max_iterations = 500
    macro_cycles = 5
    alternate_nonbonded_off_on = False
    rmsd_bonds_termination_cutoff = 0
    rmsd_angles_termination_cutoff = 0
    grms_termination_cutoff = 0
    correct_special_position_tolerance = 1.0
    riding_h = True
    move {
      bond = True
      nonbonded = True
      angle = True
      dihedral = True
      chirality = True
      planarity = True
      parallelity = True
    }
  }
-------------------------------------------------------------------------------
Processing inputs
*****************
  Monomer Library directory:
    "/home/wangwensong/phenix-dev-3409/modules/chem_data/mon_lib"
  Total number of atoms: 80
  Number of models: 1
  Model: ""
    Number of chains: 1
    Chain: "A"
      Number of atoms: 80
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 20, 80
          Classifications: {'peptide': 20}
          Link IDs: {'TRANS': 19}
  Time building chain proxies: 0.03, per 1000 atoms: 0.38
  Number of scatterers: 80
  At special positions: 0
  Unit cell: (69.499, 63.088, 83.652, 90, 90, 90)
  Space group: P 1 (No. 1)
  Number of sites at special positions: 0
  Number of scattering types: 3
    Type Number    sf(0)
     O      20      8.00
     N      20      7.00
     C      40      6.00
    sf(0) = scattering factor at diffraction angle 0.

  Number of disulfides: simple=0, symmetry=0
  Custom bonds:
    bond:
      atom 1: "ATOM      4  O   GLY A   1 .*.     O  "
      atom 2: "ATOM     17  N   GLY A   5 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.173
      distance_ideal:   2.900
      ideal - model:   -0.273
      slack:            0.000
      delta_slack:     -0.273
      sigma:            0.1000
    bond:
      atom 1: "ATOM      8  O   GLY A   2 .*.     O  "
      atom 2: "ATOM     21  N   GLY A   6 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.243
      distance_ideal:   2.900
      ideal - model:   -1.343
      slack:            0.000
      delta_slack:     -1.343
      sigma:            0.1000
    bond:
      atom 1: "ATOM     12  O   GLY A   3 .*.     O  "
      atom 2: "ATOM     25  N   GLY A   7 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   5.226
      distance_ideal:   2.900
      ideal - model:   -2.326
      slack:            0.000
      delta_slack:     -2.326
      sigma:            0.1000
    bond:
      atom 1: "ATOM     16  O   GLY A   4 .*.     O  "
      atom 2: "ATOM     29  N   GLY A   8 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.289
      distance_ideal:   2.900
      ideal - model:   -1.389
      slack:            0.000
      delta_slack:     -1.389
      sigma:            0.1000
    bond:
      atom 1: "ATOM     20  O   GLY A   5 .*.     O  "
      atom 2: "ATOM     33  N   GLY A   9 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.534
      distance_ideal:   2.900
      ideal - model:   -1.634
      slack:            0.000
      delta_slack:     -1.634
      sigma:            0.1000
    bond:
      atom 1: "ATOM     24  O   GLY A   6 .*.     O  "
      atom 2: "ATOM     37  N   GLY A  10 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.190
      distance_ideal:   2.900
      ideal - model:   -0.290
      slack:            0.000
      delta_slack:     -0.290
      sigma:            0.1000
    bond:
      atom 1: "ATOM     28  O   GLY A   7 .*.     O  "
      atom 2: "ATOM     41  N   GLY A  11 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   2.352
      distance_ideal:   2.900
      ideal - model:    0.548
      slack:            0.000
      delta_slack:      0.548
      sigma:            0.1000
    bond:
      atom 1: "ATOM     32  O   GLY A   8 .*.     O  "
      atom 2: "ATOM     45  N   GLY A  12 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.324
      distance_ideal:   2.900
      ideal - model:   -0.424
      slack:            0.000
      delta_slack:     -0.424
      sigma:            0.1000
    bond:
      atom 1: "ATOM     36  O   GLY A   9 .*.     O  "
      atom 2: "ATOM     49  N   GLY A  13 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   2.770
      distance_ideal:   2.900
      ideal - model:    0.130
      slack:            0.000
      delta_slack:      0.130
      sigma:            0.1000
    bond:
      atom 1: "ATOM     40  O   GLY A  10 .*.     O  "
      atom 2: "ATOM     53  N   GLY A  14 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   2.716
      distance_ideal:   2.900
      ideal - model:    0.184
      slack:            0.000
      delta_slack:      0.184
      sigma:            0.1000
    bond:
      atom 1: "ATOM     44  O   GLY A  11 .*.     O  "
      atom 2: "ATOM     57  N   GLY A  15 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.002
      distance_ideal:   2.900
      ideal - model:   -1.102
      slack:            0.000
      delta_slack:     -1.102
      sigma:            0.1000
    bond:
      atom 1: "ATOM     48  O   GLY A  12 .*.     O  "
      atom 2: "ATOM     61  N   GLY A  16 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.196
      distance_ideal:   2.900
      ideal - model:   -1.296
      slack:            0.000
      delta_slack:     -1.296
      sigma:            0.1000
    bond:
      atom 1: "ATOM     52  O   GLY A  13 .*.     O  "
      atom 2: "ATOM     65  N   GLY A  17 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.402
      distance_ideal:   2.900
      ideal - model:   -0.502
      slack:            0.000
      delta_slack:     -0.502
      sigma:            0.1000
    bond:
      atom 1: "ATOM     56  O   GLY A  14 .*.     O  "
      atom 2: "ATOM     69  N   GLY A  18 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.903
      distance_ideal:   2.900
      ideal - model:   -1.003
      slack:            0.000
      delta_slack:     -1.003
      sigma:            0.1000
    bond:
      atom 1: "ATOM     60  O   GLY A  15 .*.     O  "
      atom 2: "ATOM     73  N   GLY A  19 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.527
      distance_ideal:   2.900
      ideal - model:   -0.627
      slack:            0.000
      delta_slack:     -0.627
      sigma:            0.1000
    bond:
      atom 1: "ATOM     64  O   GLY A  16 .*.     O  "
      atom 2: "ATOM     77  N   GLY A  20 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.175
      distance_ideal:   2.900
      ideal - model:   -1.275
      slack:            0.000
      delta_slack:     -1.275
      sigma:            0.1000
    Total number of custom bonds: 16
  Custom angles:
    angle:
      atom 1: "ATOM      4  O   GLY A   1 .*.     O  "
      atom 2: "ATOM     17  N   GLY A   5 .*.     N  "
      atom 3: "ATOM      3  C   GLY A   1 .*.     C  "
      angle_model:   12.35
      angle_ideal:  147.15
      ideal - model:   134.80
      sigma: 5
    angle:
      atom 1: "ATOM      8  O   GLY A   2 .*.     O  "
      atom 2: "ATOM     21  N   GLY A   6 .*.     N  "
      atom 3: "ATOM      7  C   GLY A   2 .*.     C  "
      angle_model:   11.81
      angle_ideal:  147.15
      ideal - model:   135.34
      sigma: 5
    angle:
      atom 1: "ATOM     12  O   GLY A   3 .*.     O  "
      atom 2: "ATOM     25  N   GLY A   7 .*.     N  "
      atom 3: "ATOM     11  C   GLY A   3 .*.     C  "
      angle_model:   11.91
      angle_ideal:  147.15
      ideal - model:   135.24
      sigma: 5
    angle:
      atom 1: "ATOM     16  O   GLY A   4 .*.     O  "
      atom 2: "ATOM     29  N   GLY A   8 .*.     N  "
      atom 3: "ATOM     15  C   GLY A   4 .*.     C  "
      angle_model:   14.40
      angle_ideal:  147.15
      ideal - model:   132.75
      sigma: 5
    angle:
      atom 1: "ATOM     20  O   GLY A   5 .*.     O  "
      atom 2: "ATOM     33  N   GLY A   9 .*.     N  "
      atom 3: "ATOM     19  C   GLY A   5 .*.     C  "
      angle_model:   15.21
      angle_ideal:  147.15
      ideal - model:   131.94
      sigma: 5
    angle:
      atom 1: "ATOM     24  O   GLY A   6 .*.     O  "
      atom 2: "ATOM     37  N   GLY A  10 .*.     N  "
      atom 3: "ATOM     23  C   GLY A   6 .*.     C  "
      angle_model:   19.90
      angle_ideal:  147.15
      ideal - model:   127.25
      sigma: 5
    angle:
      atom 1: "ATOM     28  O   GLY A   7 .*.     O  "
      atom 2: "ATOM     41  N   GLY A  11 .*.     N  "
      atom 3: "ATOM     27  C   GLY A   7 .*.     C  "
      angle_model:   11.28
      angle_ideal:  147.15
      ideal - model:   135.87
      sigma: 5
    angle:
      atom 1: "ATOM     32  O   GLY A   8 .*.     O  "
      atom 2: "ATOM     45  N   GLY A  12 .*.     N  "
      atom 3: "ATOM     31  C   GLY A   8 .*.     C  "
      angle_model:   12.74
      angle_ideal:  147.15
      ideal - model:   134.41
      sigma: 5
    angle:
      atom 1: "ATOM     36  O   GLY A   9 .*.     O  "
      atom 2: "ATOM     49  N   GLY A  13 .*.     N  "
      atom 3: "ATOM     35  C   GLY A   9 .*.     C  "
      angle_model:   11.70
      angle_ideal:  147.15
      ideal - model:   135.45
      sigma: 5
    angle:
      atom 1: "ATOM     40  O   GLY A  10 .*.     O  "
      atom 2: "ATOM     53  N   GLY A  14 .*.     N  "
      atom 3: "ATOM     39  C   GLY A  10 .*.     C  "
      angle_model:   10.08
      angle_ideal:  147.15
      ideal - model:   137.07
      sigma: 5
    angle:
      atom 1: "ATOM     44  O   GLY A  11 .*.     O  "
      atom 2: "ATOM     57  N   GLY A  15 .*.     N  "
      atom 3: "ATOM     43  C   GLY A  11 .*.     C  "
      angle_model:   15.99
      angle_ideal:  147.15
      ideal - model:   131.16
      sigma: 5
    angle:
      atom 1: "ATOM     48  O   GLY A  12 .*.     O  "
      atom 2: "ATOM     61  N   GLY A  16 .*.     N  "
      atom 3: "ATOM     47  C   GLY A  12 .*.     C  "
      angle_model:   15.34
      angle_ideal:  147.15
      ideal - model:   131.81
      sigma: 5
    angle:
      atom 1: "ATOM     52  O   GLY A  13 .*.     O  "
      atom 2: "ATOM     65  N   GLY A  17 .*.     N  "
      atom 3: "ATOM     51  C   GLY A  13 .*.     C  "
      angle_model:    9.28
      angle_ideal:  147.15
      ideal - model:   137.87
      sigma: 5
    angle:
      atom 1: "ATOM     56  O   GLY A  14 .*.     O  "
      atom 2: "ATOM     69  N   GLY A  18 .*.     N  "
      atom 3: "ATOM     55  C   GLY A  14 .*.     C  "
      angle_model:   10.17
      angle_ideal:  147.15
      ideal - model:   136.98
      sigma: 5
    angle:
      atom 1: "ATOM     60  O   GLY A  15 .*.     O  "
      atom 2: "ATOM     73  N   GLY A  19 .*.     N  "
      atom 3: "ATOM     59  C   GLY A  15 .*.     C  "
      angle_model:   14.10
      angle_ideal:  147.15
      ideal - model:   133.05
      sigma: 5
    angle:
      atom 1: "ATOM     64  O   GLY A  16 .*.     O  "
      atom 2: "ATOM     77  N   GLY A  20 .*.     N  "
      atom 3: "ATOM     63  C   GLY A  16 .*.     C  "
      angle_model:   15.12
      angle_ideal:  147.15
      ideal - model:   132.03
      sigma: 5
    Total number of custom angles: 16
  Custom bonds:
    bond:
      atom 1: "ATOM      4  O   GLY A   1 .*.     O  "
      atom 2: "ATOM     17  N   GLY A   5 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.173
      distance_ideal:   2.900
      ideal - model:   -0.273
      slack:            0.000
      delta_slack:     -0.273
      sigma:            0.1000
    bond:
      atom 1: "ATOM      8  O   GLY A   2 .*.     O  "
      atom 2: "ATOM     21  N   GLY A   6 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.243
      distance_ideal:   2.900
      ideal - model:   -1.343
      slack:            0.000
      delta_slack:     -1.343
      sigma:            0.1000
    bond:
      atom 1: "ATOM     12  O   GLY A   3 .*.     O  "
      atom 2: "ATOM     25  N   GLY A   7 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   5.226
      distance_ideal:   2.900
      ideal - model:   -2.326
      slack:            0.000
      delta_slack:     -2.326
      sigma:            0.1000
    bond:
      atom 1: "ATOM     16  O   GLY A   4 .*.     O  "
      atom 2: "ATOM     29  N   GLY A   8 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.289
      distance_ideal:   2.900
      ideal - model:   -1.389
      slack:            0.000
      delta_slack:     -1.389
      sigma:            0.1000
    bond:
      atom 1: "ATOM     20  O   GLY A   5 .*.     O  "
      atom 2: "ATOM     33  N   GLY A   9 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.534
      distance_ideal:   2.900
      ideal - model:   -1.634
      slack:            0.000
      delta_slack:     -1.634
      sigma:            0.1000
    bond:
      atom 1: "ATOM     24  O   GLY A   6 .*.     O  "
      atom 2: "ATOM     37  N   GLY A  10 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.190
      distance_ideal:   2.900
      ideal - model:   -0.290
      slack:            0.000
      delta_slack:     -0.290
      sigma:            0.1000
    bond:
      atom 1: "ATOM     28  O   GLY A   7 .*.     O  "
      atom 2: "ATOM     41  N   GLY A  11 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   2.352
      distance_ideal:   2.900
      ideal - model:    0.548
      slack:            0.000
      delta_slack:      0.548
      sigma:            0.1000
    bond:
      atom 1: "ATOM     32  O   GLY A   8 .*.     O  "
      atom 2: "ATOM     45  N   GLY A  12 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.324
      distance_ideal:   2.900
      ideal - model:   -0.424
      slack:            0.000
      delta_slack:     -0.424
      sigma:            0.1000
    bond:
      atom 1: "ATOM     36  O   GLY A   9 .*.     O  "
      atom 2: "ATOM     49  N   GLY A  13 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   2.770
      distance_ideal:   2.900
      ideal - model:    0.130
      slack:            0.000
      delta_slack:      0.130
      sigma:            0.1000
    bond:
      atom 1: "ATOM     40  O   GLY A  10 .*.     O  "
      atom 2: "ATOM     53  N   GLY A  14 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   2.716
      distance_ideal:   2.900
      ideal - model:    0.184
      slack:            0.000
      delta_slack:      0.184
      sigma:            0.1000
    bond:
      atom 1: "ATOM     44  O   GLY A  11 .*.     O  "
      atom 2: "ATOM     57  N   GLY A  15 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.002
      distance_ideal:   2.900
      ideal - model:   -1.102
      slack:            0.000
      delta_slack:     -1.102
      sigma:            0.1000
    bond:
      atom 1: "ATOM     48  O   GLY A  12 .*.     O  "
      atom 2: "ATOM     61  N   GLY A  16 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.196
      distance_ideal:   2.900
      ideal - model:   -1.296
      slack:            0.000
      delta_slack:     -1.296
      sigma:            0.1000
    bond:
      atom 1: "ATOM     52  O   GLY A  13 .*.     O  "
      atom 2: "ATOM     65  N   GLY A  17 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.402
      distance_ideal:   2.900
      ideal - model:   -0.502
      slack:            0.000
      delta_slack:     -0.502
      sigma:            0.1000
    bond:
      atom 1: "ATOM     56  O   GLY A  14 .*.     O  "
      atom 2: "ATOM     69  N   GLY A  18 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.903
      distance_ideal:   2.900
      ideal - model:   -1.003
      slack:            0.000
      delta_slack:     -1.003
      sigma:            0.1000
    bond:
      atom 1: "ATOM     60  O   GLY A  15 .*.     O  "
      atom 2: "ATOM     73  N   GLY A  19 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   3.527
      distance_ideal:   2.900
      ideal - model:   -0.627
      slack:            0.000
      delta_slack:     -0.627
      sigma:            0.1000
    bond:
      atom 1: "ATOM     64  O   GLY A  16 .*.     O  "
      atom 2: "ATOM     77  N   GLY A  20 .*.     N  "
      symmetry operation: x,y,z
      distance_model:   4.175
      distance_ideal:   2.900
      ideal - model:   -1.275
      slack:            0.000
      delta_slack:     -1.275
      sigma:            0.1000
    Total number of custom bonds: 16
  Custom angles:
    angle:
      atom 1: "ATOM      4  O   GLY A   1 .*.     O  "
      atom 2: "ATOM     17  N   GLY A   5 .*.     N  "
      atom 3: "ATOM      3  C   GLY A   1 .*.     C  "
      angle_model:   12.35
      angle_ideal:  147.15
      ideal - model:   134.80
      sigma: 5
    angle:
      atom 1: "ATOM      8  O   GLY A   2 .*.     O  "
      atom 2: "ATOM     21  N   GLY A   6 .*.     N  "
      atom 3: "ATOM      7  C   GLY A   2 .*.     C  "
      angle_model:   11.81
      angle_ideal:  147.15
      ideal - model:   135.34
      sigma: 5
    angle:
      atom 1: "ATOM     12  O   GLY A   3 .*.     O  "
      atom 2: "ATOM     25  N   GLY A   7 .*.     N  "
      atom 3: "ATOM     11  C   GLY A   3 .*.     C  "
      angle_model:   11.91
      angle_ideal:  147.15
      ideal - model:   135.24
      sigma: 5
    angle:
      atom 1: "ATOM     16  O   GLY A   4 .*.     O  "
      atom 2: "ATOM     29  N   GLY A   8 .*.     N  "
      atom 3: "ATOM     15  C   GLY A   4 .*.     C  "
      angle_model:   14.40
      angle_ideal:  147.15
      ideal - model:   132.75
      sigma: 5
    angle:
      atom 1: "ATOM     20  O   GLY A   5 .*.     O  "
      atom 2: "ATOM     33  N   GLY A   9 .*.     N  "
      atom 3: "ATOM     19  C   GLY A   5 .*.     C  "
      angle_model:   15.21
      angle_ideal:  147.15
      ideal - model:   131.94
      sigma: 5
    angle:
      atom 1: "ATOM     24  O   GLY A   6 .*.     O  "
      atom 2: "ATOM     37  N   GLY A  10 .*.     N  "
      atom 3: "ATOM     23  C   GLY A   6 .*.     C  "
      angle_model:   19.90
      angle_ideal:  147.15
      ideal - model:   127.25
      sigma: 5
    angle:
      atom 1: "ATOM     28  O   GLY A   7 .*.     O  "
      atom 2: "ATOM     41  N   GLY A  11 .*.     N  "
      atom 3: "ATOM     27  C   GLY A   7 .*.     C  "
      angle_model:   11.28
      angle_ideal:  147.15
      ideal - model:   135.87
      sigma: 5
    angle:
      atom 1: "ATOM     32  O   GLY A   8 .*.     O  "
      atom 2: "ATOM     45  N   GLY A  12 .*.     N  "
      atom 3: "ATOM     31  C   GLY A   8 .*.     C  "
      angle_model:   12.74
      angle_ideal:  147.15
      ideal - model:   134.41
      sigma: 5
    angle:
      atom 1: "ATOM     36  O   GLY A   9 .*.     O  "
      atom 2: "ATOM     49  N   GLY A  13 .*.     N  "
      atom 3: "ATOM     35  C   GLY A   9 .*.     C  "
      angle_model:   11.70
      angle_ideal:  147.15
      ideal - model:   135.45
      sigma: 5
    angle:
      atom 1: "ATOM     40  O   GLY A  10 .*.     O  "
      atom 2: "ATOM     53  N   GLY A  14 .*.     N  "
      atom 3: "ATOM     39  C   GLY A  10 .*.     C  "
      angle_model:   10.08
      angle_ideal:  147.15
      ideal - model:   137.07
      sigma: 5
    angle:
      atom 1: "ATOM     44  O   GLY A  11 .*.     O  "
      atom 2: "ATOM     57  N   GLY A  15 .*.     N  "
      atom 3: "ATOM     43  C   GLY A  11 .*.     C  "
      angle_model:   15.99
      angle_ideal:  147.15
      ideal - model:   131.16
      sigma: 5
    angle:
      atom 1: "ATOM     48  O   GLY A  12 .*.     O  "
      atom 2: "ATOM     61  N   GLY A  16 .*.     N  "
      atom 3: "ATOM     47  C   GLY A  12 .*.     C  "
      angle_model:   15.34
      angle_ideal:  147.15
      ideal - model:   131.81
      sigma: 5
    angle:
      atom 1: "ATOM     52  O   GLY A  13 .*.     O  "
      atom 2: "ATOM     65  N   GLY A  17 .*.     N  "
      atom 3: "ATOM     51  C   GLY A  13 .*.     C  "
      angle_model:    9.28
      angle_ideal:  147.15
      ideal - model:   137.87
      sigma: 5
    angle:
      atom 1: "ATOM     56  O   GLY A  14 .*.     O  "
      atom 2: "ATOM     69  N   GLY A  18 .*.     N  "
      atom 3: "ATOM     55  C   GLY A  14 .*.     C  "
      angle_model:   10.17
      angle_ideal:  147.15
      ideal - model:   136.98
      sigma: 5
    angle:
      atom 1: "ATOM     60  O   GLY A  15 .*.     O  "
      atom 2: "ATOM     73  N   GLY A  19 .*.     N  "
      atom 3: "ATOM     59  C   GLY A  15 .*.     C  "
      angle_model:   14.10
      angle_ideal:  147.15
      ideal - model:   133.05
      sigma: 5
    angle:
      atom 1: "ATOM     64  O   GLY A  16 .*.     O  "
      atom 2: "ATOM     77  N   GLY A  20 .*.     N  "
      atom 3: "ATOM     63  C   GLY A  16 .*.     C  "
      angle_model:   15.12
      angle_ideal:  147.15
      ideal - model:   132.03
      sigma: 5
    Total number of custom angles: 16

  Automatic linking
    Parameters for automatic linking
      Linking & cutoffs
        Metal                : False - 3.50
        Amimo acid           : False - 1.90
        Carbohydrate         : True  - 1.99
        Ligands              : True  - 1.99
        Small molecules      : False - 1.98
        Amino acid - RNA/DNA : False
      
  Number of custom bonds: simple=0, symmetry=0
  Time building additional restraints: 0.01
  Conformation dependent library (CDL) restraints added in 4.8 milliseconds
  
  Adding C-beta torsion restraints...
  Number of C-beta restraints generated:  0

  Time building geometry restraints manager: 0.04 seconds

  NOTE: a complete listing of the restraints can be obtained by requesting
        output of .geo file.

  Histogram of bond lengths:
        1.20 -     1.27: 20
        1.27 -     1.34: 10
        1.34 -     1.41: 9
        1.41 -     1.48: 18
        1.48 -     1.55: 22
  Bond restraints: 79
  Sorted by residual:
  bond pdb=" CA  GLY A   3 "
       pdb=" C   GLY A   3 "
    ideal  model  delta    sigma   weight residual
    1.514  1.455  0.059 1.41e-02 5.03e+03 1.76e+01
  bond pdb=" C   GLY A   5 "
       pdb=" O   GLY A   5 "
    ideal  model  delta    sigma   weight residual
    1.235  1.198  0.037 1.35e-02 5.49e+03 7.40e+00
  bond pdb=" N   GLY A   8 "
       pdb=" CA  GLY A   8 "
    ideal  model  delta    sigma   weight residual
    1.453  1.421  0.032 1.24e-02 6.50e+03 6.67e+00
  bond pdb=" CA  GLY A  15 "
       pdb=" C   GLY A  15 "
    ideal  model  delta    sigma   weight residual
    1.510  1.546 -0.036 1.67e-02 3.59e+03 4.64e+00
  bond pdb=" N   GLY A   1 "
       pdb=" CA  GLY A   1 "
    ideal  model  delta    sigma   weight residual
    1.451  1.484 -0.033 1.60e-02 3.91e+03 4.13e+00
  ... (remaining 74 not shown)

  Histogram of bond angle deviations from ideal:
      108.05 -   111.36: 8
      111.36 -   114.67: 12
      114.67 -   117.98: 13
      117.98 -   121.28: 43
      121.28 -   124.59: 21
  Bond angle restraints: 97
  Sorted by residual:
  angle pdb=" N   GLY A  11 "
        pdb=" CA  GLY A  11 "
        pdb=" C   GLY A  11 "
      ideal   model   delta    sigma   weight residual
     114.67  111.35    3.32 1.10e+00 8.26e-01 9.10e+00
  angle pdb=" O   GLY A  11 "
        pdb=" C   GLY A  11 "
        pdb=" N   GLY A  12 "
      ideal   model   delta    sigma   weight residual
     122.84  119.46    3.38 1.16e+00 7.43e-01 8.51e+00
  angle pdb=" CA  GLY A   9 "
        pdb=" C   GLY A   9 "
        pdb=" N   GLY A  10 "
      ideal   model   delta    sigma   weight residual
     117.80  114.61    3.19 1.11e+00 8.12e-01 8.28e+00
  angle pdb=" CA  GLY A   8 "
        pdb=" C   GLY A   8 "
        pdb=" O   GLY A   8 "
      ideal   model   delta    sigma   weight residual
     120.66  117.70    2.96 1.06e+00 8.90e-01 7.78e+00
  angle pdb=" CA  GLY A   9 "
        pdb=" C   GLY A   9 "
        pdb=" O   GLY A   9 "
      ideal   model   delta    sigma   weight residual
     119.50  122.45   -2.95 1.18e+00 7.18e-01 6.26e+00
  ... (remaining 92 not shown)

  Histogram of dihedral angle deviations from ideal:
        0.90 -     5.64: 8
        5.64 -    10.38: 3
       10.38 -    15.12: 6
       15.12 -    19.87: 1
       19.87 -    24.61: 1
  Dihedral angle restraints: 19
    sinusoidal: 0
      harmonic: 19
  Sorted by residual:
  dihedral pdb=" CA  GLY A  16 "
           pdb=" C   GLY A  16 "
           pdb=" N   GLY A  17 "
           pdb=" CA  GLY A  17 "
      ideal   model   delta  harmonic     sigma   weight residual
     180.00  155.39   24.61     0      5.00e+00 4.00e-02 2.42e+01
  dihedral pdb=" CA  GLY A   6 "
           pdb=" C   GLY A   6 "
           pdb=" N   GLY A   7 "
           pdb=" CA  GLY A   7 "
      ideal   model   delta  harmonic     sigma   weight residual
     180.00  163.74   16.26     0      5.00e+00 4.00e-02 1.06e+01
  dihedral pdb=" CA  GLY A   8 "
           pdb=" C   GLY A   8 "
           pdb=" N   GLY A   9 "
           pdb=" CA  GLY A   9 "
      ideal   model   delta  harmonic     sigma   weight residual
     180.00  165.40   14.60     0      5.00e+00 4.00e-02 8.53e+00
  ... (remaining 16 not shown)

  Chirality restraints: 0

  Planarity restraints: 19
  Sorted by residual:
                                 delta    sigma   weight rms_deltas residual
  plane pdb=" CA  GLY A   5 "    0.013 2.00e-02 2.50e+03   2.64e-02 6.97e+00
        pdb=" C   GLY A   5 "   -0.046 2.00e-02 2.50e+03
        pdb=" O   GLY A   5 "    0.017 2.00e-02 2.50e+03
        pdb=" N   GLY A   6 "    0.015 2.00e-02 2.50e+03
                                 delta    sigma   weight rms_deltas residual
  plane pdb=" CA  GLY A  17 "   -0.012 2.00e-02 2.50e+03   2.33e-02 5.45e+00
        pdb=" C   GLY A  17 "    0.040 2.00e-02 2.50e+03
        pdb=" O   GLY A  17 "   -0.015 2.00e-02 2.50e+03
        pdb=" N   GLY A  18 "   -0.013 2.00e-02 2.50e+03
                                 delta    sigma   weight rms_deltas residual
  plane pdb=" CA  GLY A   9 "   -0.011 2.00e-02 2.50e+03   2.21e-02 4.87e+00
        pdb=" C   GLY A   9 "    0.038 2.00e-02 2.50e+03
        pdb=" O   GLY A   9 "   -0.014 2.00e-02 2.50e+03
        pdb=" N   GLY A  10 "   -0.013 2.00e-02 2.50e+03
  ... (remaining 16 not shown)

  Histogram of nonbonded interaction distances:
        2.57 -     3.03: 41
        3.03 -     3.48: 45
        3.48 -     3.93: 83
        3.93 -     4.38: 110
        4.38 -     4.84: 144
  Nonbonded interactions: 423
  Sorted by model distance:
  nonbonded pdb=" N   GLY A   5 "
            pdb=" N   GLY A   6 "
     model   vdw
     2.573 2.560
  nonbonded pdb=" O   GLY A  11 "
            pdb=" CA  GLY A  12 "
     model   vdw
     2.632 2.752
  nonbonded pdb=" N   GLY A   1 "
            pdb=" N   GLY A   2 "
     model   vdw
     2.641 2.560
  nonbonded pdb=" N   GLY A   6 "
            pdb=" N   GLY A   7 "
     model   vdw
     2.660 2.560
  nonbonded pdb=" O   GLY A  15 "
            pdb=" CA  GLY A  16 "
     model   vdw
     2.661 2.752
  ... (remaining 418 not shown)

  NOTE: a complete listing of the restraints can be obtained by requesting
        output of .geo file.
-------------------------------------------------------------------------------
Atom selection
**************
  selected 80 atoms out of total 80
-------------------------------------------------------------------------------
Geometry Restraints
*******************
-------------------------------------------------------------------------------
Setup riding H
**************
-------------------------------------------------------------------------------
Minimization
************
    target: 13727.3
      bond_residual_sum (n=95): 1988.79
      nonbonded_residual_sum (n=423): 0.946699
      angle_residual_sum (n=113): 11627.2
      dihedral_residual_sum (n=19): 82.7736
      chirality_residual_sum (n=0): 0
      planarity_residual_sum (n=19): 27.5874
      parallelity_residual_sum (n=0): 0
  macro-cycle: 0
    target: 9469.68
      bond_residual_sum (n=95): 364.038
      nonbonded_residual_sum (n=842): 32.0929
      angle_residual_sum (n=113): 8903.48
      dihedral_residual_sum (n=19): 154.543
      chirality_residual_sum (n=0): 0
      planarity_residual_sum (n=19): 15.5322
      parallelity_residual_sum (n=0): 0
  macro-cycle: 1
    target: 9418.27
      bond_residual_sum (n=95): 382.747
      nonbonded_residual_sum (n=842): 30.4864
      angle_residual_sum (n=113): 8868.28
      dihedral_residual_sum (n=19): 123.359
      chirality_residual_sum (n=0): 0
      planarity_residual_sum (n=19): 13.3945
      parallelity_residual_sum (n=0): 0
  macro-cycle: 2
    target: 9418.27
      bond_residual_sum (n=95): 382.778
      nonbonded_residual_sum (n=841): 30.4843
      angle_residual_sum (n=113): 8868.24
      dihedral_residual_sum (n=19): 123.372
      chirality_residual_sum (n=0): 0
      planarity_residual_sum (n=19): 13.3954
      parallelity_residual_sum (n=0): 0
  macro-cycle: 3
    target: 9418.27
      bond_residual_sum (n=95): 382.778
      nonbonded_residual_sum (n=841): 30.4843
      angle_residual_sum (n=113): 8868.24
      dihedral_residual_sum (n=19): 123.372
      chirality_residual_sum (n=0): 0
      planarity_residual_sum (n=19): 13.3954
      parallelity_residual_sum (n=0): 0
  macro-cycle: 4
    target: 9418.27
      bond_residual_sum (n=95): 382.778
      nonbonded_residual_sum (n=841): 30.4843
      angle_residual_sum (n=113): 8868.24
      dihedral_residual_sum (n=19): 123.372
      chirality_residual_sum (n=0): 0
      planarity_residual_sum (n=19): 13.3954
      parallelity_residual_sum (n=0): 0
-------------------------------------------------------------------------------
Write PDB file
**************
  output file name: m-helix_shaken_minimized.pdb
min,max,mean shift from start:  1.066  6.280  3.366
min,max,mean shift from start:  1.066  6.280  3.366
-------------------------------------------------------------------------------
Write GEO file
**************
-------------------------------------------------------------------------------
Model statistics
****************

Geometry Restraints Library: GeoStd + Monomer Library + CDL v1.2
Deviations from Ideal Values.
  Bond      :  0.011   0.046     79
  Angle     :  2.014   5.060     97
  Chirality :  0.000   0.000      0
  Planarity :  0.008   0.011     19
  Dihedral  : 12.741  15.854     19
  Min Nonbonded Distance : 2.557

Molprobity Statistics.
  All-atom Clashscore : 830.99
  Ramachandran Plot:
    Outliers :  0.00 %
    Allowed  :  0.00 %
    Favored  : 100.00 %
  Rotamer Outliers :  0.00 %
  Cbeta Deviations :  0.00 %
  Peptide Plane:
    Cis-proline     : 0.00 %
    Cis-general     : 0.00 %
    Twisted Proline : 0.00 %
    Twisted General : 0.00 %
-------------------------------------------------------------------------------
Detailed timing
***************
  Initialization, inputs:   0.110
  Processing inputs:        0.450
  Atom selection:           0.000
  Geometry Restraints:      0.000
  Setup riding H:           0.000
  Minimization:             0.420
  Write PDB file:           0.000
  Write GEO file:           0.050
  Model statistics:         0.720
  Sum of individual times: 1.750
Overall runtime: 1.750   
