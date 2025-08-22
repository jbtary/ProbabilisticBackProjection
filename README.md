# Probabilistic Back-projection algorithm 

In this repository are included the Matlab codes needed to carry out the probabilistic back-projection of an event. These codes were originally developed for the study Li and Gudmundsson, 2020, and extended for the study Tary et al., 2025. 

To locate an event in 2D with a constant velocity, follow loc_whales_2d_oneevent.m
To locate an event in 3D with a constant velocity, follow loc_whales_3d_oneevent.m
To locate an event in 3D with a travel-time tables, follow loc_whales_tt_oneevent.m

To make synthetic tests in 2D, 3D or using travel-time tables, follow synthetic_2d.m, synthetic_3d.m, and synthetic_tt.m.

To use the above-mentioned codes, you have to build a mapping operator that can be prepared and saved in advance using the codes build_mapping_operator_2d.m, build_mapping_operator_3d.m, or build_mapping_operator_tt.m.

While the function to transform cross-correlation envelopes to probabilites can be constructed (see make_distri_stat.m in prob_back_loc folder), it can be calculated from the data as well (see distri_stat.m and gen_tran_func_rat_func.m).

Three examples of events (whale calls) are included in the Events_mat folder. They can be plotted using the seeVAevent.m.

No travel-time tables and velocity model examples are included here due to their large sizes. Refer to the documentation of NonLinLoc for indications as to how to create them (Vel2Grid and Grid2Times only).

References to be cited in publications using these codes:

Li, K. L., & Gudmundsson, O. (2020). A probabilistic tremor location method. Geophysical Research Letters, 47(4), e2019GL085538.

Tary, J. B., Poveda, S. F., Li, K. L., Peirce, C., Hobbs, R. W., & Vargas, C. A. (2025). Detection and localization of Bryde’s whale calls using machine learning and probabilistic back-projection. The Journal of the Acoustical Society of America, 158(2), 1386-1397.
