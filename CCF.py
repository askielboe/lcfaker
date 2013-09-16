import numpy as np
import matplotlib.pyplot as plt

class CCF():
    def __init__(self, time, ccf):
        self.time = np.array(time)
        self.ccf = np.array(ccf)
        self.threshold = 0.8
        self.multimodalThreshold = 0.8

    def getLag(self, threshold=-1):
        """
        Returns the median time-lag above some threshold
        """

        if threshold == -1:
            threshold = self.threshold

        try:
            mask = self.ccf >= threshold*max(self.ccf)
            ccfCut = self.ccf[mask]
            timesCut = self.time[mask]
        except ValueError:
            return -9999

        # Figure out if we have several peaks above threshold
        try:
            mask = self.ccf >= self.multimodalThreshold*max(self.ccf)
            timesAboveCut = self.time[mask]
            timesAboveCutSub = timesAboveCut[1:].astype(float) - timesAboveCut[:-1].astype(float)
            if np.any(timesAboveCutSub > 5.):
                return -9999
        except ValueError:
            return -9999

        # If we don't have a well defined maximum in the CFF return -9999
        # try:
        #     mask = self.ccf >= 0.7*max(self.ccf)
        #     timesAboveCut = self.time[mask]
        #     print "max-min = ", max(timesAboveCut) - min(timesAboveCut)
        #     if max(timesAboveCut) - min(timesAboveCut) > 20:
        #         return -9999
        # except ValueError:
        #     return -9999

        # Compute center using 2nd degree polynomial fit to CFF above 0.4
        #fit = np.polyfit(timesCut, ccfCut, 2)
        #center = -fit[1]/(2.*fit[0])
        #print center

        # Compute center using accumulated CCF values above cut
        # sumCCF = np.sum(ccfCut)
        # sumRunning = 0.
        # for i in range(len(ccfCut)):
        #     sumRunning += ccfCut[i]
        #     if sumRunning >= sumCCF / 2.:
        #         center = (timesCut[i] + timesCut[i-1]) / 2.
        #         break
        #print 'center = ', center

        # Compute centroid
        centroid = 0.0
        ccfCut /= np.sum(ccfCut)
        for i in range(len(ccfCut)):
            centroid += ccfCut[i] * timesCut[i]
        #print 'centroid = ', centroid

        try:
            return centroid
        except UnboundLocalError:
            return -9999
            print 'WARNING: Could not find CCF center! - Returning 0.0.'
            #raise ValueError('ERROR: Could not find CCF center!')

    def plot(self):
        plt.figure()
        plt.title("Cross Correlation Function")
        plt.xlabel("time-lag")
        plt.ylabel("CCF")
        plt.grid(True)
        plt.plot(self.time, self.ccf, label="CCF")
        plt.plot([self.time[0], self.time[-1]],
            [self.threshold*max(self.ccf), self.threshold*max(self.ccf)],'--y',
            label='Threshold')
        plt.plot([self.time[0], self.time[-1]],
            [self.multimodalThreshold*max(self.ccf), self.multimodalThreshold*max(self.ccf)],'--r',
            label='Multimodal Threshold')
        plt.legend()
        plt.show()

    def __repr__(self):
        return len(np.concatenate((self.time, self.ccf)))
