import numpy as np
import xarray as xr
import os, os.path
import glob


home_dir="/home/Feiyu.Lu/Documents/SPEAR_ECDA/"
vftmp_dir="/vftmp/Feiyu.Lu/"
arch_dir="/archive/Feiyu.Lu/"

class Exp:
    def __init__(self, dict, read_ens=False):
        """
        Initialize the experiment class with the given parameters.
        """
        self.user = dict.get('user', 'Feiyu.Lu')
        self.subfolder = dict.get('subfolder', '')
        
        if 'exp_name' in dict:
            self.exp_name = dict['exp_name']
            self.ocean_res = self.exp_name.split('_')[2]
            self.atmos_res = self.exp_name.split('_')[1]
            self.scenario = '_'.join(self.exp_name.split('_')[3:-1])
            self.exp_suffix = self.exp_name.split('_')[-1]
        else:
            self.exp_suffix = dict.get('exp_suffix', 'H64')
            self.ocean_res = dict.get('ocean_res', 'o1')
            self.atmos_res = dict.get('atmos_res', 'c96')
            self.scenario = dict.get('scenario', 'ECDA')
            self.exp_name='SPEAR_{}_{}_{}_{}'.format(self.atmos_res,self.ocean_res,self.scenario,self.exp_suffix)
            
        self.exp_dir='/archive/{}/SPEAR/{}/{}'.format(self.user,self.subfolder,self.exp_name)
        print(self.exp_dir)

        if (os.path.exists('{}/ensemble'.format(self.exp_dir)) or 
            os.path.exists('{}/pp_ensemble'.format(self.exp_dir)) or 
            os.path.exists('{}/ens_01'.format(self.exp_dir)) or 
            os.path.exists('{}/pp_ens_01'.format(self.exp_dir))):
            self.ensemble = True
        else:
            self.ensemble = False

        if 'years' in dict:
            self.years = dict['years']
        else:
            if os.path.exists('{}/ensemble'.format(self.exp_dir)):
                ensemble_done = sorted(glob.glob('{}/ensemble/*.done'.format(self.exp_dir)))
                self.years = [int(done.split('/')[-1].split('.')[0][0:4]) for done in ensemble_done]

            elif os.path.exists('{}/ens_01'.format(self.exp_dir)) or os.path.exists('{}/pp_ens_01'.format(self.exp_dir)):
                ens_01_done = sorted(glob.glob('{}/ens_01/*.done'.format(self.exp_dir)))
                self.years = [int(done.split('/')[-1].split('.')[0][0:4]) for done in ens_01_done]

            else:
                history_files = sorted(glob.glob('{}/history/*.nc.tar'.format(self.exp_dir)))
                self.years = [int(file.split('/')[-1].split('.')[0][0:4]) for file in history_files]
                
        if len(self.years) == self.years[-1] - self.years[0] + 1:
            self.year_start = self.years[0]
            self.year_end = self.years[-1]
        else:
            raise ValueError('Years are not continuous.')
            
        self.output_dir='{}{}/'.format(home_dir,self.exp_name)
        if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)

        self.read_ens = read_ens
        self.ens_member = []
        if self.ensemble:
            self.pp_ensemble = os.path.exists('{}/pp_ensemble'.format(self.exp_dir))
            if self.pp_ensemble:
                self.ens_member.append(Member('{}/pp_ensemble'.format(self.exp_dir)))

            self.pp_ens = os.path.exists('{}/pp_ens_01'.format(self.exp_dir))
            self.ensemble_size = max(len(sorted(glob.glob('{}/ens_??'.format(self.exp_dir)))),
                                     len(sorted(glob.glob('{}/pp_ens_??'.format(self.exp_dir)))))
            if self.pp_ens and self.read_ens:
                for i in range(self.ensemble_size):
                    self.ens_member.append(Member('{}/pp_ens_{:02d}'.format(self.exp_dir,i+1)))
        else:
            self.ensemble_size = 1
            self.pp = os.path.exists('{}/pp'.format(self.exp_dir))
            if self.pp:
                self.ens_member.append(Member('{}/pp'.format(self.exp_dir)))

    def get_years(self):
        return self.years
    
    def get_ensemble_size(self):
        return self.ensemble_size
        
    def get_component_names(self, type=[], ensemble_number=0):
        if type:
            comps = [comp for comp in self.ens_member[ensemble_number].component_names if comp.split('_')[0] in type]
        else:
            comps = self.ens_member[ensemble_number].component_names
        return comps
    
    def get_ts_freq(self, component_name, ensemble_number=0):
        return self.ens_member[ensemble_number].components[component_name].ts_freq
    
    def get_ts_variable_names(self, component_name, freq, ensemble_number=0):
        return self.ens_member[ensemble_number].components[component_name].ts_variables[freq]
    
    def get_av_freq(self, component_name, ensemble_number=0):
        return self.ens_member[ensemble_number].components[component_name].av_freq
    
    def get_ts_files(self, component_name, freq, variable_name, years=[], ensemble_number=0):
        files = self.ens_member[ensemble_number].components[component_name].ts_files[freq][variable_name]
        if len(years) > 0:
            files = [file for file in files if int(file.split('/')[-1].split('.')[-3][0:4]) in years]
        return files
    
    def get_ts_files_vftmp(self, component_name, freq, variable_name, years=[], ensemble_number=0):
        files = self.ens_member[ensemble_number].components[component_name].ts_files[freq][variable_name]
        if len(years) > 0:
            files = [file for file in files if int(file.split('/')[-1].split('.')[-3][0:4]) in years]
        for file in files:
            vftmp_file = file.replace(arch_dir,vftmp_dir)
            if os.path.exists(vftmp_file):
                files[files.index(file)] = vftmp_file
            else:
                raise ValueError('File {} does not exist in vftmp'.format(vftmp_file))
        return files
    
    def get_ts_ds(self, component_name, freq, variable_name, years=[], ensemble_number=0):
        try:
            ds = xr.open_mfdataset(self.get_ts_files_vftmp(component_name, freq, variable_name, years, ensemble_number))
            print('Using vftmp files')
        except:
            ds = xr.open_mfdataset(self.get_ts_files(component_name, freq, variable_name, years, ensemble_number))
            print('Using archive files')
        return ds
    
    def get_ts_files_ens(self, component_name, freq, variable_name, years=[], ensemble_numbers=[]):
        if ensemble_numbers == []:
            ensemble_numbers = range(1,self.ensemble_size+1)
        files = []
        for ensemble_number in ensemble_numbers:
            ens_files = self.ens_member[ensemble_number].components[component_name].ts_files[freq][variable_name]
            if len(years) > 0:
                ens_files = [file for file in ens_files if int(file.split('/')[-1].split('.')[-3][0:4]) in years]
            files.append(ens_files)
        return files
    
    def get_ts_ds_ens(self, component_name, freq, variable_name, years=[], ensemble_numbers=[]):
        if ensemble_numbers == []:
            ensemble_numbers = range(1,self.ensemble_size+1)
        try:
            ds = xr.open_mfdataset(self.get_ts_files_ens_vftmp(component_name, freq, variable_name, years, ensemble_numbers),
                               concat_dim=['ens','time'], combine='nested')
            print('Using vftmp files')
        except:
            ds = xr.open_mfdataset(self.get_ts_files_ens(component_name, freq, variable_name, years, ensemble_numbers),
                                concat_dim=['ens','time'], combine='nested')
            print('Using archive files')
        ds["ens"] = ("ens", ensemble_numbers)
        return ds
    
    def get_av_files(self, component_name, freq, years=[], ensemble_number=0):
        files = self.ens_member[ensemble_number].components[component_name].av_files[freq]
        if len(years) > 0:
            files = [file for file in files if int(file.split('/')[-1].split('.')[-3][0:4]) in years]
        return files
    
    def get_av_ds(self, component_name, freq, years=[], ensemble_number=0):
        ds = xr.open_mfdataset(self.get_av_files(component_name, freq, years, ensemble_number))
        return ds
    
    def get_ts_files_ens_vftmp(self, component_name, freq, variable_name, years=[], ensemble_numbers=[]):
        if ensemble_numbers == []:
            ensemble_numbers = range(1,self.ensemble_size+1)

        files = []
        for ensemble_number in ensemble_numbers:
            ens_files = self.ens_member[ensemble_number].components[component_name].ts_files[freq][variable_name]
            if len(years) > 0:
                ens_files = [file for file in ens_files if int(file.split('/')[-1].split('.')[-3][0:4]) in years]
            for file in ens_files:
                vftmp_file = file.replace(arch_dir,vftmp_dir)
                if os.path.exists(vftmp_file):
                    ens_files[ens_files.index(file)] = vftmp_file
                else:
                    raise ValueError('File {} does not exist in vftmp'.format(vftmp_file))
                    
            files.append(ens_files)

        return files
    
    def get_ts_ds_ens_vftmp(self, component_name, freq, variable_name, years=[], ensemble_numbers=[]):
        if ensemble_numbers == []:
            ensemble_numbers = range(1,self.ensemble_size+1)
        ds = xr.open_mfdataset(self.get_ts_files_ens_vftmp(component_name, freq, variable_name, years, ensemble_numbers),
                               concat_dim=['ens','time'], combine='nested')
        ds["ens"] = ("ens", ensemble_numbers)
        return ds
    
class Member:
    def __init__(self, member_dir):
        """
        Initialize the member class with the given parameters.
        """
        self.member_dir = member_dir
        component_names = os.listdir(self.member_dir)
        self.component_names = []
        self.components = {}
        for component_name in component_names:
            if component_name.split('_')[0] in ['atmos','ocean','ice','land']:
                self.component_names.append(component_name)
                self.components[component_name]=Component('{}/{}'.format(self.member_dir,component_name),component_name)
                if not self.components[component_name].ts_freq:
                    del self.components[component_name]
            
class Component:
    def __init__(self, component_dir, component_name):
        """
        Initialize the component class with the given parameters.
        """
        self.component_dir = component_dir
        self.component_name = component_name
        try:
            self.static_file = glob.glob('{}/{}*static.nc'.format(self.component_dir,component_name))[0]
        except:
            self.static_file = []
        
        self.variable_names = []
        self.ts_freq = []
        self.ts_files = {}
        self.ts_dirs = {}
        self.ts_variables = {}

        if os.path.exists('{}/ts'.format(self.component_dir)):
            self.ts_dir = '{}/ts'.format(self.component_dir)

            for freq, length in zip(['daily','monthly','monthly','annual'],[1,1,10,10]):
                if os.path.exists('{}/ts/{}/{}yr'.format(self.component_dir,freq,length)):
                    ts_name='{}_{}yr'.format(freq,length)
                    self.ts_freq.append(ts_name)
                    self.ts_dirs[ts_name] = '{}/ts/{}/{}yr'.format(self.component_dir,freq,length)
                    self.ts_variables[ts_name] = list(set([file.split('/')[-1].split('.')[-2] 
                                                           for file in os.listdir(self.ts_dirs[ts_name])]))
                    self.ts_files[ts_name] = {}
                    for var in self.ts_variables[ts_name]:
                        self.ts_files[ts_name][var] = \
                            sorted(glob.glob('{}/*.{}.nc'.format(self.ts_dirs[ts_name],var)))

        if os.path.exists('{}/av'.format(self.component_dir)):
            self.av_dir = '{}/av'.format(self.component_dir)
            self.av_freq = []
            self.av_files = {}
            for av_file_dir in os.listdir(self.av_dir):
                if av_file_dir.split('_')[0] in ['annual','monthly']:
                    files = sorted(glob.glob('{}/{}/*.nc'.format(self.av_dir,av_file_dir)))
                    if len(files) > 0:
                        self.av_freq.append(av_file_dir)
                        self.av_files[av_file_dir] = \
                            sorted(glob.glob('{}/{}/*.nc'.format(self.av_dir,av_file_dir)))
            