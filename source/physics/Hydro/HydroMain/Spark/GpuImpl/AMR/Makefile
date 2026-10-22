Hydro += hy_gpu_amr_types.o hy_gpu_amr_kernels.o hy_gpu_amr_advance.o hy_gpu_amr_adapter.o
hy_gpu_amr_kernels.o: hy_gpu_amr_types.o
hy_gpu_amr_advance.o: hy_gpu_amr_kernels.o
hy_gpu_amr_adapter.o: hy_gpu_amr_advance.o Grid_interface.o Grid_tile.o Grid_iterator.o Hydro_data.o Eos_data.o Eos_interface.o hy_rk_interface.o Driver_interface.o RuntimeParameters_interface.o
Hydro.o: hy_gpu_amr_adapter.o
