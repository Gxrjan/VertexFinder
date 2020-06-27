using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

/*
    This class is responsible for reading a polytope from a file
*/
class PolytopeReader
{
    /// <summary>
    /// Reads a file and returns polytope in cdd H-format
    /// </summary>
    /// <param name="src">File path</param>
    /// <returns>Polytope</returns>
    internal static Polytope read_from_file_cdd(string src)
    {
        Polytope p = null;
        try
        {
            List<Inequality> inequalitiesList = new List<Inequality>();
            StreamReader sr = new StreamReader(src);
            string line;
            while ((line = sr.ReadLine()) != null)
            {
                string[] parts = Regex.Split(line.Trim(), "[ ]+");
                double[] coefficients = new double[parts.Length];
                coefficients[coefficients.Length - 1] = Convert.ToDouble(parts[0]);
                for (int i = 1; i < parts.Length; i++)
                {
                    coefficients[i - 1] = -1 * Convert.ToDouble(parts[i]);
                }
                inequalitiesList.Add(new Inequality(coefficients));
            }
            p = new Polytope(inequalitiesList.ToArray());


        }
        catch (Exception E)
        {
            Console.WriteLine("File not found");
        }
        return p;
    }
}
