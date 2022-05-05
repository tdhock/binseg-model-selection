from ruptures.costs import NotEnoughPoints
from ruptures.base import BaseCost
import numpy as np
from numpy.linalg import norm

class CostL2CumSum(BaseCost):

    r"""
    Least squared deviation.

    For two indexes $a<b$, the cost function on the sub-signal
    $y_{a..b} = [y_a, y_{a+1},\dots,y_{b-1}]$, is equal to

    $$
     c(y_{a..b}) = \sum_{t=a}^{b-1} \|y_t - \bar{y}_{a..b}\|^2
    $$

    where $\bar{y}_{a..b}$ is the empirical mean of $y_{a..b}$.

    For efficiency, the cost function is re-written with the cumulative sums of
    the signals $y$ and $[\|y_1\|^2, \|y_2\|^2,\dots,\|y_T\|^2]$.
    """

    model = "l2_cumsum"

    def __init__(self):
        """Initialize the object."""
        self.signal = None
        self.signal_cumsum = None
        self.signal_norm_cumsum = None
        self.min_size = 1

    def fit(self, signal) -> "CostL2":
        """Set parameters of the instance.

        Args:
            signal (array): array of shape (n_samples,) or (n_samples, n_features)

        Returns:
            self
        """
        if signal.ndim == 1:
            self.signal = signal.reshape(-1, 1)
        else:
            self.signal = signal

        self.signal_cumsum = np.cumsum(self.signal, axis=0)
        self.signal_norm_cumsum = np.cumsum(norm(self.signal, axis=1) ** 2, axis=0)

        return self

    def error(self, start, end) -> float:
        """Return the approximation cost on the segment [start:end].

        Args:
            start (int): start of the segment
            end (int): end of the segment

        Returns:
            segment cost (float): the value of the cost on the segment [start:end]

        Raises:
            NotEnoughPoints: when the segment is too short (less than `min_size`
            samples).
        """
        if end - start < self.min_size:
            raise NotEnoughPoints

        if start == 0:
            res = self.signal_norm_cumsum[end - 1] - norm(
                self.signal_cumsum[end - 1]
            ) ** 2 / (end - start)
            return res
        res = (
            self.signal_norm_cumsum[end - 1]
            - self.signal_norm_cumsum[start - 1]
            - norm(self.signal_cumsum[end - 1] - self.signal_cumsum[start - 1]) ** 2
            / (end - start)
        )
        return res
