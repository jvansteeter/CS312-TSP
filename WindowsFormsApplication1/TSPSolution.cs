using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TSP;

namespace WindowsFormsApplication1
{
	class TSPSolution
	{
		/// <summary>
		/// we use the representation [cityB,cityA,cityC] 
		/// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
		/// and the edge from cityC to cityB is the final edge in the path.  
		/// You are, of course, free to use a different representation if it would be more convenient or efficient 
		/// for your data structure(s) and search algorithm. 
		/// </summary>
		public ArrayList
			Route;

		/// <summary>
		/// constructor
		/// </summary>
		/// <param name="iroute">a (hopefully) valid tour</param>
		public TSPSolution(ArrayList iroute)
		{
			Route = new ArrayList(iroute);
		}

		/// <summary>
		/// Compute the cost of the current route.  
		/// Note: This does not check that the route is complete.
		/// It assumes that the route passes from the last city back to the first city. 
		/// </summary>
		/// <returns></returns>
		public double costOfRoute()
		{
			// go through each edge in the route and add up the cost. 
			int x;
			City here;
			double cost = 0D;

			for (x = 0; x < Route.Count - 1; x++)
			{
				here = Route[x] as City;
				cost += here.costToGetTo(Route[x + 1] as City);
			}

			// go from the last city to the first. 
			here = Route[Route.Count - 1] as City;
			cost += here.costToGetTo(Route[0] as City);
			return cost;
		}
	}
}
