import math 
import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap

## 2D Ising model 
# initial state E0 generated randomly
up = 1
down = -1 
B = 0.5
J = 100 
T = 2
size = 100 
MC_steps = 200
initial = np.random.choice([-1,1], size= (size,size))

# thermodynamics proporties of Ising plotting 

def Ham(state): # obtain the internal energy of the system 
    mag = - B * np.sum(state)
    cor  =  0 
    for i in range(size): 
        for j in range(size):
            nn_sum = state[(i+1)% size, j]+ state[(i-1)%size , j]+state[i,(j+1)%size]+ state[i,(j-1)%size]
            cor -= J*  nn_sum
    return mag+cor

def entropy (state): # S = kB ln W
    # unsure what to use so i use the stat defintion'kB ln W ' on a 2D lattice
    # in this case kB is set to be 1 by default 
    # may also implied sterling's approx: 
    # ln x!/(x-m)!
    n_up = (size**2 +abs(np.sum(state)))/2
    tot = size**2
    return tot*((n_up/tot)*math.log(n_up/tot) +(1-n_up/tot)*math.log(1-n_up/tot))

def helm (state): # A= U-TS
    return Ham(state)- T*entropy(state)
def gibbs (state): # G = A +BM
    return helm(state) + B*np.sum(state)
def mag(state): 
    return np.sum(state)/(size**2)
    
# H = -B * s_i - sum_ J* si* sj
# if -1 to + 1 flip, mag term changes : -2 B; correlation change -2J
# if 1 to -1, mag change +2B: correlation change +2J 
# for spin in state[i,j]

# code for performing MC and the visulisation of the lattice state
fig,(ax1,ax2)  = plt.subplots(1,2)
cmap = ListedColormap(['white', 'black'])
MC_step_list = range(MC_steps)
#basic settings for MC steps and plotting
m_list = [mag(initial)]
# some arrays to store the thermo variables
# change mag to other physical varibales defined above to visualize other quantities 
im = ax1.imshow(initial, cmap=cmap, interpolation='nearest')
pl = ax2.scatter(MC_step_list[0], mag(initial), color = 'g', label ='average magnetisation', s = 1)
ax2.set(xlim= [0,MC_steps], ylim = [-1, 1])
ax2.legend()
data = initial 

def MC_step (state): 
    for _ in range( size* size): 
    # select a random spin
        [i,j]= np.random.randint(0, size, 2 )
    # calculate the energy change if a flip has been done 
        nn_sum = state[(i+1)% size, j]+ state[(i-1)%size , j]+state[i,(j+1)%size]+ state[i,(j-1)%size]
        d_E =  B * 2 * state[i,j] +2* J * state[i,j]* nn_sum
        if d_E< 0 or np.exp(-d_E/T) > 1/size**2: 
            state[i,j]*= -1 
    return state

def update(frame): 
    # plot the Magnetisation over steps 
    x = MC_step_list[:frame]
    data = MC_step(initial)
    m_list.append(mag(data))
    # change mag to other physical varibales defined above to visualize other quantities 
    y = m_list[:frame]
   
    pl.set_offsets( np.stack([x,y]).T)
    im.set_array(data)
    return (pl,im)

# Create the animation
ani = FuncAnimation(fig=fig, func = update, frames=MC_steps, interval=100)
#ani.save('simulation_with_real_time_m.gif')

plt.show()
    