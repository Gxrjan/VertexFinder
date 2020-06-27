using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using System.Threading;
using Google.OrTools.LinearSolver;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Double.Solvers;
using MathNet.Numerics.LinearAlgebra.Solvers;

class Program
{

    public static void Main(string[] args)
    {
        Polytope p;
        String path = "";
        try
        {
            int k = -1;
            switch (args.Length)
            {
                case 0:
                    Console.Write("Enter file path: ");
                    path = Console.ReadLine();
                    break;
                case 1:
                    path = args[0];
                    start(path);
                    break;
                case 3:
                    path = args[0];
                    string option = args[1];
                    int num = int.Parse(args[2]);
                    if (option == "-e")
                        start_k(path, num);
                    else if (option == "-v")
                        start_v(path, num);
                    else Console.WriteLine("Wrong usage"); Environment.Exit(0);
                    k = int.Parse(args[1]);
                    break;
                default:
                    Console.WriteLine("Wrong usage");
                    Environment.Exit(0);
                    break;
            }
        }
        catch (Exception E)
        {
            Console.WriteLine(E.Message);
        }

    }


    /// <summary>
    /// Starts the process of exploring the face lattice till the user stops the program
    /// </summary>
    /// <param name="path">Path to the file with a polytope</param>
    private static void start(string path)
    {

        try
        {
            Polytope p = PolytopeReader.read_from_file_cdd(path);
            int iteration = 1;
            while (true)
            {
                p.explore(iteration, p.vertices);
                iteration++;
            }
        }
        catch (Exception E)
        {
            Console.WriteLine(E.Message);
        }

    }


    /// <summary>
    /// Start the process of exploring the face lattice untill the program finds a certain number of vertices
    /// </summary>
    /// <param name="path">Path to the file with a polytope</param>
    /// <param name="num">The program will explore untill it finds this number of vertices</param>
    private static void start_v(string path, int num)
    {
        try
        {
            Polytope p = PolytopeReader.read_from_file_cdd(path);
            int iteration = 1;
            while (p.vertices.Count < num)
            {
                p.explore(iteration, p.vertices);
                iteration++;
            }
        } catch (IOException e)
        {
            Console.WriteLine(e.Message);
        }
    }



    /// <summary>
    /// The program will terminate exploring if it can't find another vertex in n^k iterations, 
    /// where n is the current number of vertices
    /// and k is the number from the command line arguments
    /// </summary>
    /// <param name="path">Path to the file with a polytope</param>
    /// <param name="num">k</param>
    private static void start_k(string path, int num)
    {
        try
        {
            Polytope p = PolytopeReader.read_from_file_cdd(path);
            int iteration = 1;
            int number_of_vertices = 0;
            int last_discovery_iteration = 0;
            while (true)
            {
                p.explore(iteration, p.vertices);
                if (p.vertices.Count > number_of_vertices)
                {
                    number_of_vertices = p.vertices.Count;
                    last_discovery_iteration = iteration;
                }
                if (iteration - last_discovery_iteration > (number_of_vertices ^ num))
                {
                    break;
                }
                iteration++;
            }

        } catch (IOException e)
        {
            Console.WriteLine(e.Message);
        }
    }
}