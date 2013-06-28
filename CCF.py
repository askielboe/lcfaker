import numpy as np
import matplotlib.pyplot as plt

class CCF():
    def __init__(self, time, ccf):
        self.time = np.array(time)
        self.ccf = np.array(ccf)

    def getLag(self, threshold=0.8):
        """
        Returns the median time-lag above some threshold
        """
        mask = self.ccf >= threshold*max(self.ccf)
        ccfCut = self.ccf[mask]
        timesCut = self.time[mask]

        # Compute center using 2nd degree polynomial fit to CFF above 0.4
        #fit = np.polyfit(timesCut, ccfCut, 2)
        #center = -fit[1]/(2.*fit[0])
        #print center

        # Compute center using accumulated CCF values above cut of 0.2
        sumCCF = np.sum(ccfCut)
        sumRunning = 0.
        for i in range(len(ccfCut)):
            sumRunning += ccfCut[i]
            if sumRunning >= sumCCF / 2.:
                center = (timesCut[i] + timesCut[i-1]) / 2.
                break

        try:
            return center
        except UnboundLocalError:
            return 0.0
            print 'WARNING: Could not find CCF center! - Returning 0.0.'
            #raise ValueError('ERROR: Could not find CCF center!')

    def plot(self):
        plt.figure()
        plt.title("Cross Correlation Function")
        plt.xlabel("time-lag")
        plt.ylabel("CCF")
        plt.grid(True)
        plt.plot(self.time, self.ccf)
        plt.show()

    def __repr__(self):
        return len(np.concatenate((self.time, self.ccf)))
