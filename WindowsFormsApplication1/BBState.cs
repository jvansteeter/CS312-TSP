using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WindowsFormsApplication1
{
    class BBState : IComparable<BBState>
    {
        public double[,] rcm;
		public double lowerBound;
		public List<int> route;
		public static double arnellConstant;

        public BBState(double[,] rcm, double lowerBound, List<int> route)
        {
            this.rcm = rcm;
            this.lowerBound = lowerBound;
            this.route = route;
        }

		public int CompareTo(BBState other)
		{
			// to determine priority score = depth^2 / lowerBound
			// highest score has highest priority
			double myScore = Math.Pow(route.Count, (/*(route.Count + 1) / */Math.Sqrt(route.Count))) / lowerBound;
			double otherScore = Math.Pow(other.route.Count, (/*(other.route.Count + 1) / */Math.Sqrt(other.route.Count))) / other.lowerBound;
			if (myScore > otherScore)
				return -1;
			else
				return 1;
		}

		/*public int CompareTo(BBState other)
		{
			// to determine priority score = depth^2 / lowerBound
			// highest score has highest priority
			double myScore = route.Count * route.Count / lowerBound;
			double otherScore = other.route.Count * other.route.Count / other.lowerBound;
			if (myScore > otherScore)
				return -1;
			else
				return 1;
		}*/
    }
}
