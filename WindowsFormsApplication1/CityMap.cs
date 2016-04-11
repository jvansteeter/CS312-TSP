using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TSP;

namespace WindowsFormsApplication1
{
	class CityMap
	{
		int[][] map;

		// a map of the closest cities to a given city
		public CityMap(City[] Cities, int numNearby)
		{
			map = new int[Cities.Length][];

			for(int y = 0; y < Cities.Length; y++)
			{
				int[] closestCities = new int[numNearby];
				List<double> distances = new List<double>();
				for(int x = 0; x < Cities.Length; x++)
				{
					distances.Add(Cities[y].costToGetTo(Cities[x]));
				}
				for (int i = 0; i < closestCities.Length; i++)
				{
					closestCities[i] = distances.IndexOf(distances.Min());
					distances[distances.IndexOf(distances.Min())] = Double.PositiveInfinity;
				}
				map[y] = closestCities;
			}
		}

		public int[] getClosestCities(int city)
		{
			return map[city];
		}
	}
}
