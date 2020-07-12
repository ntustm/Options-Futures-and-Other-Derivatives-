import numpy as np
import matplotlib.pyplot as plt

fo = open('ko_barrier.log', 'w')
def log_nd(nd, name='null'):
    global fo
    log = ''
    for l in nd:
        log += str(l.tolist())
        log += '\n'
    fo.write('%s = [\n%s]\n' % (name, log))

#paras

nperiod = 8
nday_aperiod = 5
S0 = 48970
para = dict(
        S0 = S0,
        K = S0*0.97,
        H = S0*1.0186,
        T_days = nperiod*nday_aperiod,
        obs_days_lst = (np.arange(nperiod)+1)*nday_aperiod, #np.linspace(1,nperiod,nperiod)*nday_aperiod,
        r = 0.03,
        q = 0.03,
        sigma = 0.13,
        num_tests = 50000,
        call_put = 'call',
        up_down = 'up'      
)
#paras

def McGbmQ(S0,r,sigma,T,num_tests,nStep):
    # monte-carlo, Geometric Brownian motion
    # rng(10)
    W = np.random.randn(num_tests,nStep)
    #print(W)
    #plt.subplot(311)
    #plt.plot(W)
    h = T/nStep # dt  resolution 1/245
    #fo.write("dt=%.f\n" %  h)
    #log_nd(W)
    dlogS = (r-sigma*sigma/2)*h + sigma*np.sqrt(h)*W
    #log_nd(np.cumsum(dlogS,1))
    #plt.subplot(312)
    #plt.plot(dlogS)
    S = S0 * np.exp(np.cumsum(dlogS,1))
    fo.write("")
    #plt.subplot(313)
    #for l in S:
    #    l=np.insert(l, 0, S0)
    #    plt.plot(l)
    #plt.show()
    return S


def knockoutbarrier(S0,K,H,T_days,obs_days_lst,r,q,sigma,num_tests,call_put='call',up_down='up'):
    #print(obs_days_lst)
    T = T_days/243.0
    nStep = T_days
    dt = T/nStep
    S = McGbmQ(S0,r-q,sigma,T,num_tests,nStep)
    price_lst = [0]*num_tests
    for i in range(0,num_tests):
        knock_out_triggered = 0
        for j in range(len(obs_days_lst)):
            if up_down == 'up':
                if S[i][obs_days_lst[j]-1] > H:
                    price_lst[i] = 0
                    knock_out_triggered  =1
                    break
            else:
                if S[i][obs_days_lst[j]-1] < H:
                    price_lst[i] = 0
                    knock_out_triggered  =1
                    break
        if knock_out_triggered == 0:
            if call_put == 'call':
                price_lst[i] = np.exp(-r*T)*max(S[i][-1] - K,0)
            else:
                price_lst[i] = np.exp(-r*T)*max(K - S[i][-1],0)
    return np.mean(price_lst)

def calc_ko_price():
    from scipy.optimize import fsolve
    def ff(H):
        S0 = 49000
        K = 0.97
        ntest = 50000
        sigma = 0.13
        return sum([
            (knockoutbarrier(S0,S0*K,S0*H,i*5,(np.arange(i)+1)*5,0.03,0.03,sigma,ntest,'call','up')
             -1.8* knockoutbarrier(S0,S0*K,S0*H,i*5,(np.arange(i)+1)*5,0.03,0.03,sigma,ntest,'put','up'))
        for i in range(1,9)])
    return fsolve(ff,0.99) 


if __name__ == "__main__":
    print(knockoutbarrier(**para))
    print(calc_ko_price())
