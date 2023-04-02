#-----------------------------------------------------------------------------#
# Author: Mohamed Watfa
# URL: https://github.com/medwatt/
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Changes: Barak Hoffer
# URL: https://github.com/barakhoffer/
#-----------------------------------------------------------------------------#
import numpy as np
import ngspyce as ns
import os

def range_to_arr(r):
    start, stop, step = r
    return np.arange(start, stop + step, step)


class LookupTableGenerator:
    def __init__(
        self,
        vgs=(0, 1.8, 0.05),
        vds=(0, 1, 0.5),
        vsb=(0, 1, 0.5),
        width=1.0,
        lengths=(0.15, 0.3, 0.5),
        temp=27,
        model_names={'nmos' : 'sky130_fd_pr__nfet_01v8', 'pmos' : 'sky130_fd_pr__pfet_01v8'},
        description="gmid lookup table",
    ):
        self.vgs = np.array(vgs)
        self.vds = np.array(vds)
        self.vsb = np.array(vsb)
        self.width = width
        self.lengths = np.array(lengths)
        self.temp = temp
        self.model_names = model_names
        self.description = description
        self.pdk_path = os.environ['PDKPATH']
        self.parameter_names = [
            "id",
            "vth",
            "gm",
            "gmbs",
            "gds",
            "vdsat",
            "cgg",
            "cgs",
            "cbg",
            "cgd",
            "cdd",
        ]
        self.lookup_table = {}

    def __initalize(self, mos):
        self.n_lengths = len(self.lengths)
        self.n_vsb = round((self.vsb[1] - self.vsb[0]) / self.vsb[2]) + 1
        self.n_vds = round((self.vds[1] - self.vds[0]) / self.vds[2]) + 1
        self.n_vgs = round((self.vgs[1] - self.vgs[0]) / self.vgs[2]) + 1
        self.lookup_table[mos] = {}
        for p in self.parameter_names:
            self.lookup_table[mos][p] = np.zeros(
                shape=(self.n_lengths, self.n_vsb, self.n_vds, self.n_vgs)
            )


    def __load_circuit(self, mos, default_length):
        ns.cmd('reset')
        ns.cmd('destroy all')
        # ToDo make more generic includes
        ns.circ([".param mc_mm_switch=0",
                 ".param mc_pr_switch=0",
                f".include {self.pdk_path}/libs.tech/ngspice/corners/tt.spice",
                f".include {self.pdk_path}/libs.tech/ngspice/r+c/res_typical__cap_typical.spice",
                f".include {self.pdk_path}/libs.tech/ngspice/r+c/res_typical__cap_typical__lin.spice",
                f".include {self.pdk_path}/libs.tech/ngspice/corners/tt/specialized_cells.spice",
                 "vgs ng gnd 0",
                 "vbs nb gnd dc 0",
                 "vds nd gnd 0",
                f".param length={default_length}",
                f"x{mos} nd ng nb gnd {self.model_names[mos]} l={{length}} w={self.width} nf=1 ad='int((nf+1)/2) * W/nf * 0.29' as='int((nf+2)/2) * W/nf * 0.29' pd='2*int((nf+1)/2) * (W/nf + 0.29)' ps='2*int((nf+2)/2) * (W/nf + 0.29)' nrd='0.29 / W' nrs='0.29 / W' sa=0 sb=0 sd=0 mult=1 m=1"])
        ns.cmd('run')

    def __simulate_circuit(self, mos, length, vsb):
        r = 1 if mos == 'nmos' else -1
        ns.alterparams(length=length)
        ns.alter('vbs', dc=-vsb * r)
        [ns.save(f"@m.x{mos}.m{self.model_names[mos]}[{p}]") for p in self.parameter_names]
        dc_sim = ['vds'] + list(self.vds * r)
        dc_sim += ['vgs'] + list(self.vgs * r)
        analysis = ns.dc(*dc_sim)
        return analysis

    def __save_parameters(self, analysis, mos, length, vsb):
        for p in self.parameter_names:
            res = analysis[f"@m.x{mos}.m{self.model_names[mos]}[{p}]"]
            self.lookup_table[mos][p][length][vsb] = res

    def __generate_loopkup_table(self, mos):
        self.__initalize(mos)
        self.__load_circuit(mos, self.lengths[0])
        for idx, length in enumerate(self.lengths):
            for idy, vsb in enumerate(
                np.linspace(self.vsb[0], self.vsb[1], self.n_vsb)
            ):
                analysis = self.__simulate_circuit(mos, length, vsb)
                self.__save_parameters(analysis, mos, idx, idy)

    def __save_to_dictionary(self):
        self.lookup_table["nmos"]["vgs"] = range_to_arr(self.vgs)
        self.lookup_table["nmos"]["vds"] = range_to_arr(self.vds)
        self.lookup_table["nmos"]["vsb"] = range_to_arr(self.vsb)
        self.lookup_table["nmos"]["w"] = self.width
        self.lookup_table["nmos"]["l"] = self.lengths
        self.lookup_table["nmos"]["parameter_names"] = self.parameter_names
        self.lookup_table["nmos"]["model_name"] = self.model_names["nmos"]
        self.lookup_table["pmos"]["vgs"] = -range_to_arr(self.vgs)
        self.lookup_table["pmos"]["vds"] = -range_to_arr(self.vds)
        self.lookup_table["pmos"]["vsb"] = -range_to_arr(self.vsb)
        self.lookup_table["pmos"]["w"] = self.width
        self.lookup_table["pmos"]["l"] = self.lengths
        self.lookup_table["pmos"]["parameter_names"] = self.parameter_names
        self.lookup_table["pmos"]["model_name"] = self.model_names["pmos"]
        self.lookup_table["description"] = self.description

    def build(self, filepath):
        # NMOS
        print("Generating lookup table for NMOS")
        self.__generate_loopkup_table('nmos')
        # PMOS
        print("Generating lookup table for PMOS")
        self.__generate_loopkup_table('pmos')
        # save results to dictionary and to file
        print("Saving to file")
        self.__save_to_dictionary()
        np.save(filepath, self.lookup_table, allow_pickle=True)
        return self.lookup_table


if __name__ == "__main__":

    # example
    obj = LookupTableGenerator (
        vgs=(0, 1, 0.01),
        vds=(0, 1, 0.05),
        vsb=(0, 1, 0.1),
        width=10e-6,
        lengths=[50e-9, 100e-9, 200e-9, 400e-9, 800e-9, 1.6e-6, 3.2e-6, 6.4e-6],
        models_path="./models",
        model_names={
            "nmos": "NMOS_VTH",
            "pmos": "PMOS_VTH"},
        description="freepdk 45nm"
        )
    obj.build("./freepdk45_loopup_table.npy")
