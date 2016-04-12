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
				List<KeyValuePair<int, double>> distances = new List<KeyValuePair<int, double>>();
				for(int x = 0; x < Cities.Length; x++)
				{
					distances.Add(new KeyValuePair<int,double>(x, Cities[y].costToGetTo(Cities[x])));
				}
				for(int i = 0; i < distances.Count; i++)
				{
					for(int j = 0; j < distances.Count - 1; j++)
					{
						if(distances[j].Value > distances[j + 1].Value)
						{
							KeyValuePair<int, double> temp = distances[j];
							distances[j] = distances[j + 1];
							distances[j + 1] = temp;
						}
					}
				}
				for (int i = 0; i < distances.Count; i++)
				{
					closestCities[i] = distances[i].Key;
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
