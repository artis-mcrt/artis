import numpy as np
import matplotlib.pyplot as plt
import sys
from os import listdir
from os.path import isfile, join

# define some parameters
plt.rcParams['xtick.labelsize']=4
plt.rcParams['ytick.labelsize']=4
plt.rcParams['axes.labelsize']=4
plt.rcParams['axes.titlesize']=4

# class for data we use for timing
class timing_data:
    # declare datatype
    def __init__(self):
        self.rank = np.array([],dtype=np.int32)
        self.thread = np.array([],dtype=np.int32)
        self.duration = np.array([],dtype=np.float32)
        self.start = np.array([],dtype=np.float32)
        self.end = np.array([],dtype=np.float32)
        self.index = np.array([],dtype=np.float32)
        self.function = []
        self.nthreads = 0
        self.nranks = 0
        self.size = 0
        self.sim_start = 0

    # Read in data from text files
    def read_outputs(self,input_file):
        # Temporary arrays to store the data
        temp_rank = []
        temp_thread = []
        temp_duration = []
        temp_start = []
        temp_end = []

        # Read a file
        with open(input_file) as f:
            for line in f:
                if "thread_timing" in line:
                    s = line.split()
                    temp_rank.append(int(s[1]))
                    temp_thread.append(int(s[2]))
                    temp_duration.append(float(s[3]))
                    temp_start.append(float(s[4]))
                    temp_end.append(float(s[5]))
                    self.function.append(s[6])

        # Convert to numpy arrays
        self.rank     = np.append(self.rank,np.asarray(temp_rank, dtype=np.int64))
        self.thread   = np.append(self.thread,np.asarray(temp_thread, dtype=np.int64))
        self.duration = np.append(self.duration,np.asarray(temp_duration, dtype=np.float64))
        self.start    = np.append(self.start,np.asarray(temp_start, dtype=np.float64))
        self.end      = np.append(self.end,np.asarray(temp_end, dtype=np.float64))
        if len(self.rank) > 0:
            self.nthreads = np.max(self.thread)+1
            self.nranks = np.max(self.rank)+1
            self.length = len(self.rank)
            self.index = self.rank*self.nthreads + self.thread

# Plot the data for a given set of files on fig
def plot_timing(files,fig,nplots,plot_num,max_x,name):
    # read the data
    data = timing_data()
    for i in files:
        data.read_outputs(i)
    
    data.sim_start = np.min(data.start)
    data.start -= data.sim_start
    data.end -= data.sim_start
    
    # Add subplot to fig
    ax = fig.add_subplot(nplots,1,plot_num)
    
    # Plot it.
    for i in range(data.length):
        c = 'b'
        if data.function[i] == 'update_packets':
            c = 'r'
        ax.barh(data.index[i],data.duration[i], left=data.start[i], height=0.95, color=c)
    
    # Messy formatting
    ax.set_xlim([0,max_x])
    labels = ['{:d} / {:d}'.format(data.rank[i], data.thread[i]) for i in range(data.length)]
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.set_yticks(data.index)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Duration (ms)")
    ax.set_ylabel("rank/thread")
    ax.set_title(name)
    fig.tight_layout()

# Find what the maximum time of each simulation is.
def get_max_x_range(dirs):
    max_x_range = 0.
    for dir in dirs:
        files = [dir+f for f in listdir(dir) if isfile(join(dir, f))]
        data = timing_data()
        for file in files:
            data.read_outputs(file)

        data.sim_start = np.min(data.start)
        data.end -= data.sim_start
        max_x_range = np.max([max_x_range,np.max(data.end)])

    return max_x_range


# Put everything together

# Create figure we're going to use
fig = plt.figure()

# Get the directories where the data is stored
nvars = len(sys.argv)
if nvars % 2 == 0 or nvars == 1:
    sys.exit("Wrong number of arguments to script. Please specify list of directories containing output files followed by list of titles.")
else:
    dirs = sys.argv[1:int((nvars-1)/2+1)]
    names = sys.argv[int((nvars-1)/2+1):]

# Find the max x range of the data
max_x_range = get_max_x_range(dirs)

# Plot the data from each directory
for i in range(len(dirs)):
    files = [dirs[i]+f for f in listdir(dirs[i]) if isfile(join(dirs[i], f))]
    plot_timing(files,fig,len(dirs),i+1, max_x_range, names[i])
    print("finished plotting "+str(i))

fig.savefig("timing.png", dpi=400)

