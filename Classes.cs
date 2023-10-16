using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Rhino;
using Rhino.Geometry;


using KangarooSolver;   //Remove if not used
using KangarooSolver.Goals;

namespace MouldLib
{
     
    ///-------SURFACE AND PIN-BED CLASSES-------///

    ///Actuator Pin-Bed Class
    public class PinBed 
    {
        //Properties currently public, make private
        public int uDiv;
        public int vDiv;
        public double uSize;
        public double vSize;
        public double maxPinHeight;
        public double lengthFactor;
        public double minCurvature;

        private Mesh bed;
        private Node[][] vertexNodes;

        //Constructors
        public PinBed(int uDiv, int vDiv, double uSize, double vSize, double maxHeight, double lengthFactor, double minCurvature)
        {
            this.uDiv = uDiv;
            this.vDiv = vDiv;
            this.uSize = uSize;
            this.vSize = vSize;
            this.maxPinHeight = maxHeight;
            this.lengthFactor = lengthFactor;
            this.minCurvature = minCurvature;

            this.bed = new Mesh();

            vertexNodes = new Node[uDiv + 1][];
            for (int i = 0; i <= uDiv; i++)
            {
                vertexNodes[i] = new Node[vDiv + 1];
            }

            CreateBed();
            SetNodeProperties();
        }

        public PinBed(Node[][] meshNodes, double maxHeight, double lengthFactor, double minCurvature)  //Creating PinBed from Node[][]
        {
            this.vertexNodes = meshNodes;

            this.uDiv = vertexNodes.Length - 1;
            this.vDiv = vertexNodes[0].Length - 1;

            Point3d m_pt = vertexNodes[0][0].GetPosition();
            Point3d test_pt = vertexNodes[1][1].GetPosition();

            this.uSize = Math.Abs(test_pt.X - m_pt.X);
            this.vSize = Math.Abs(test_pt.Y - m_pt.Y);

            this.maxPinHeight = maxHeight;
            this.lengthFactor = lengthFactor;
            this.minCurvature = minCurvature;

            CreateBedFromNodes();
            SetNodeProperties();
        }

        public PinBed(Mesh mesh, int uNum, int vNum, double maxPinHeight, double lengthFactor, double minCurvature)
        {
            this.uDiv = uNum - 1;
            this.vDiv = vNum - 1;
            this.maxPinHeight = maxPinHeight;
            this.lengthFactor = lengthFactor;
            this.minCurvature = minCurvature;
            this.bed = new Mesh();

            List<Point3d> vertexList = new List<Point3d>();
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                vertexList.Add(mesh.Vertices.Point3dAt(i));
            }

            vertexNodes = new Node[uNum][];
            for (int i = 0; i < uNum; i++)
            {
                vertexNodes[i] = new Node[vNum];
            }

            for (int i = 0; i < uNum; i++)
            {
                for (int j = 0; j < vNum; j++)
                {
                    int index = i * vNum + j;
                    vertexNodes[i][j] = new Node(vertexList[index]);
                }
            }

            Point3d m_pos = vertexNodes[0][0].GetPosition();
            Point3d next = vertexNodes[1][1].GetPosition();
            this.uSize = Math.Abs(next.X - m_pos.X);
            this.vSize = Math.Abs(next.Y - m_pos.Y);

            CreateBedFromNodes();
            SetNodeProperties();

        }

        public PinBed(List<Point3d> points, int uNum, int vNum, double maxPinHeight, double lengthFactor, double minCurvature)
        {
            this.uDiv = uNum - 1;
            this.vDiv = vNum - 1;
            this.maxPinHeight = maxPinHeight;
            this.lengthFactor = lengthFactor;
            this.minCurvature = minCurvature;
            this.bed = new Mesh();

            vertexNodes = new Node[uNum][];
            for (int i = 0; i < uNum; i++)
            {
                vertexNodes[i] = new Node[vNum];
            }

            for (int i = 0; i < uNum; i++)
            {
                for (int j = 0; j < vNum; j++)
                {
                    int index = i * vNum + j;
                    vertexNodes[i][j] = new Node(points[index]);
                }
            }

            Point3d m_pos = vertexNodes[0][0].GetPosition();
            Point3d next = vertexNodes[1][1].GetPosition();
            this.uSize = Math.Abs(next.X - m_pos.X);
            this.vSize = Math.Abs(next.Y - m_pos.Y);

            CreateBedFromNodes();
            SetNodeProperties();
        }
        //Private Methods
        private void CreateBed()  //creating mesh where vertices are also Nodes.
        {
            bed = new Mesh();
            for (int i = 0; i <= uDiv; i++)
            {
                for (int j = 0; j <= vDiv; j++)
                {
                    Point3d vert = new Point3d(i * uSize, j * vSize, 0);
                    vertexNodes[i][j] = new Node();
                    vertexNodes[i][j].SetPosition(vert);

                    bed.Vertices.Add(vertexNodes[i][j].GetPosition());
                }
            }

            for (int i = 0; i < uDiv; i++)
            {
                for (int j = 0; j < vDiv; j++)
                {
                    int a = (i * (vDiv + 1)) + j;
                    int b = (i * (vDiv + 1)) + j + 1;
                    int c = ((i + 1) * (vDiv + 1)) + j + 1;
                    int d = ((i + 1) * (vDiv + 1)) + j;

                    bed.Faces.AddFace(a, b, c, d);
                }
            }
            bed.Normals.ComputeNormals();
            bed.Compact();
        }

        private void CreateBedFromNodes()  //when preset Nodes are assigned
        {
            bed = new Mesh();
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    bed.Vertices.Add(vertexNodes[i][j].GetPosition());
                }
            }
            for (int i = 0; i < vertexNodes.Length - 1; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length - 1; j++)
                {
                    int a = (i * (vertexNodes[i].Length)) + j;
                    int b = (i * (vertexNodes[i].Length)) + j + 1;
                    int c = ((i + 1) * (vertexNodes[i].Length)) + j + 1;
                    int d = ((i + 1) * (vertexNodes[i].Length)) + j;

                    bed.Faces.AddFace(a, b, c, d);
                }
            }
            bed.Normals.ComputeNormals();
            bed.Compact();
        }

        private void SetNodeProperties()   //Similar method used in SurfaceNodes class
        {
            Rhino.Geometry.Collections.MeshVertexNormalList normals = bed.Normals;

            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    int index = i * vertexNodes[i].Length + j;  //index of current vertex
                    Vector3d normal = normals[index];   //normal at current point, can change to Vector(0,0,1)
                    Point3d pos = vertexNodes[i][j].GetPosition();  //position of current Node

                    int neighbourCount = 0;
                    double angleSum = 0.0;
                    double distSum = 0.0;

                    if (i == 0)
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = vertexNodes[i + 1][j].GetPosition();
                        Vector3d vec = neigh_pos - pos;
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec));

                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        vertexNodes[i][j].SetNeighPos(vec);
                    }

                    else if (i == vertexNodes.Length - 1)
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = vertexNodes[i - 1][j].GetPosition();
                        Vector3d vec = neigh_pos - pos;
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec));

                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        vertexNodes[i][j].SetNeighPos(vec);
                    }

                    else
                    {
                        neighbourCount += 2;
                        Point3d neigh_pos1 = vertexNodes[i + 1][j].GetPosition();
                        Point3d neigh_pos2 = vertexNodes[i - 1][j].GetPosition();

                        Vector3d vec1 = neigh_pos1 - pos;
                        Vector3d vec2 = neigh_pos2 - pos;
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec1));
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec2));

                        angleSum += Vector3d.VectorAngle(normal, vec1);
                        angleSum += Vector3d.VectorAngle(normal, vec2);
                        distSum += vec1.Length;
                        distSum += vec2.Length;

                        vertexNodes[i][j].SetNeighPos(vec1);
                        vertexNodes[i][j].SetNeighPos(vec2);

                        //Set U_Curvature
                        Plane plane = new Plane(pos, normal);
                        if (Math.Abs(plane.DistanceTo(neigh_pos1)) >= 0.001 || Math.Abs(plane.DistanceTo(neigh_pos2)) >= 0.001)
                        {
                            Circle cir = new Circle(neigh_pos2, pos, neigh_pos1);                            
                            vertexNodes[i][j].SetUCurvature(cir.Radius);
                        }
                    }

                    if (j == 0)  //Edge Points @ v = 0
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = vertexNodes[i][j + 1].GetPosition();
                        Vector3d vec = neigh_pos - pos;
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec));

                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        vertexNodes[i][j].SetNeighPos(vec);
                    }
                    else if (j == vertexNodes[i].Length - 1) //Edge Points @ v = 1
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = vertexNodes[i][j - 1].GetPosition();
                        Vector3d vec = neigh_pos - pos;
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec));
                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        vertexNodes[i][j].SetNeighPos(vec);
                    }

                    else  //middle nodes in v direction
                    {
                        neighbourCount += 2;
                        Point3d neigh_pos1 = vertexNodes[i][j + 1].GetPosition();
                        Point3d neigh_pos2 = vertexNodes[i][j - 1].GetPosition();

                        Vector3d vec1 = neigh_pos1 - pos;
                        Vector3d vec2 = neigh_pos2 - pos;
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec1));
                        vertexNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec2));

                        angleSum += Vector3d.VectorAngle(normal, vec1);
                        angleSum += Vector3d.VectorAngle(normal, vec2);
                        distSum += vec1.Length;
                        distSum += vec2.Length;

                        vertexNodes[i][j].SetNeighPos(vec1);
                        vertexNodes[i][j].SetNeighPos(vec2);

                        //Set V_Curvature
                        Plane plane = new Plane(pos, normal);
                        if (Math.Abs(plane.DistanceTo(neigh_pos1)) >= 0.001 || Math.Abs(plane.DistanceTo(neigh_pos2)) >= 0.001)
                        {
                            Circle cir = new Circle(neigh_pos2, pos, neigh_pos1);                            
                            vertexNodes[i][j].SetVCurvature(cir.Radius);
                        }
                    }

                    vertexNodes[i][j].SetAngleSum(angleSum);
                    vertexNodes[i][j].SetDistanceSum(distSum);
                    vertexNodes[i][j].SetNeighbourCount(neighbourCount);
                    vertexNodes[i][j].Set_UIndex(i);
                    vertexNodes[i][j].Set_VIndex(j);
                }
            }
        }

        //Public Methods

        public void Set_NodePanelIndices(int panel_index)
        {
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    vertexNodes[i][j].Set_PanelIndex(panel_index);
                }
            }
        }
        public Mesh GetPinBed()
        {
            return bed;
        }
        public Point3d[][] GetVertexPositions()
        {
            Point3d[][] outPoints = new Point3d[vertexNodes.Length][];
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                outPoints[i] = new Point3d[vertexNodes[i].Length];
            }

            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    outPoints[i][j] = vertexNodes[i][j].GetPosition();
                }
            }
            return outPoints;
        }
        public List<Point3d> GetPoint3dList()
        {
            List<Point3d> pts = new List<Point3d>();
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    pts.Add(vertexNodes[i][j].GetPosition());
                }
            }
            return pts;
        }
        public List<double> Get_UCurvatures()
        {
            List<double> vals = new List<double>();
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    vals.Add(vertexNodes[i][j].GetUCurvature());
                }
            }
            return vals;
        }
        public List<double> Get_VCurvatures()
        {
            List<double> vals = new List<double>();
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    vals.Add(vertexNodes[i][j].GetVCurvature());
                }
            }
            return vals;
        }
        public Node[][] GetVertexNodes()
        {
            return vertexNodes;
        }
        public void GetProperties(out int uNum, out int vNum, out double uSize, out double vSize, out double maxPinHeight, out double lengthFactor)
        {
            uNum = this.vertexNodes.Length;
            vNum = this.vertexNodes[0].Length;
            uSize = this.uSize;
            vSize = this.vSize;
            maxPinHeight = this.maxPinHeight;
            lengthFactor = this.lengthFactor;
        }
        public void MovePinBed(Vector3d motion)
        {
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    vertexNodes[i][j].MoveNode(motion);  //move all Nodes
                }
            }
            CreateBedFromNodes();
        }
        public void ApplyForce_ToPin(Vector3d force, int i_pinIndex, int j_pinIndex)
        {
            vertexNodes[i_pinIndex][j_pinIndex].AddForce(force);
        }
        public void AdjustPinBed()
        {
            for (int i = 0; i < vertexNodes.Length; i++)
            {
                for (int j = 0; j < vertexNodes[i].Length; j++)
                {
                    vertexNodes[i][j].AdjustNode();
                }
            }
            CreateBedFromNodes();
        }
        public void AdjustPinBed_ByHeight(double[][] pin_heights)
        {
            for (int i = 0; i < pin_heights.Length; i++)
            {
                for (int j = 0; j < pin_heights[i].Length; j++)
                {
                    Point3d temp = vertexNodes[i][j].GetPosition();
                    vertexNodes[i][j].SetPosition(new Point3d(temp.X, temp.Y, pin_heights[i][j]));
                }
            }
            CreateBedFromNodes();
        }
        public PinBed Duplicate()
        {
            return new PinBed(this.uDiv, this.vDiv, this.uSize, this.vSize, this.maxPinHeight, this.lengthFactor, this.minCurvature);
        }
    }

    ///Surface Nodes Class
    public class SurfaceNodes
    {
        private NurbsSurface srf;
        private int uDivs;
        private int vDivs;

        private Node[][] srfNodes;
        private double uSize;
        private double vSize;

        private Mesh mesh;

        public SurfaceNodes(Surface srf, int uDivs, int vDivs)
        {
            this.srf = srf.ToNurbsSurface();
            this.uDivs = uDivs;
            this.vDivs = vDivs;
            this.uSize = 1.0 / uDivs;
            this.vSize = 1.0 / vDivs;

            Interval interval = new Interval(0, 1);
            this.srf.SetDomain(0, interval);
            this.srf.SetDomain(1, interval);

            this.srfNodes = new Node[uDivs + 1][];
            for (int i = 0; i <= uDivs; i++)
            {
                srfNodes[i] = new Node[vDivs + 1];
            }

            SetDivisions();
            SetNodeProperties();
            CreateMesh();
        }

        //Private methods
        private void SetDivisions()
        {
            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    Node temp = new Node();
                    temp.SetPosition(srf.PointAt(i * uSize, j * vSize));
                    srfNodes[i][j] = temp;
                }
            }
        }
        private void SetNodeProperties()
        {
            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    int neighbourCount = 0;
                    double angleSum = 0.0;
                    double distSum = 0.0;

                    Vector3d normal = srf.NormalAt(i * uSize, j * vSize);  //Surface Normal at node position
                    //normal.Reverse(); //Change made after seeing results.
                    Point3d m_pos = srfNodes[i][j].GetPosition(); //Position of node

                    if (i == 0) //Edge Points @ u = 0
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = srfNodes[i + 1][j].GetPosition(); //No node in u-1 position
                        Vector3d vec = neigh_pos - m_pos;
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec)); //add vec to Node property list

                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        srfNodes[i][j].SetNeighPos(vec);
                    }
                    else if (i == srfNodes.Length - 1) //edge Points @ u = 1
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = srfNodes[i - 1][j].GetPosition(); //No node in u+1 position
                        Vector3d vec = neigh_pos - m_pos;
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec));

                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        srfNodes[i][j].SetNeighPos(vec);
                    }

                    else //middle nodes in u direction
                    {
                        neighbourCount += 2;
                        Point3d neigh_pos1 = srfNodes[i + 1][j].GetPosition();
                        Point3d neigh_pos2 = srfNodes[i - 1][j].GetPosition();

                        Vector3d vec1 = neigh_pos1 - m_pos;
                        Vector3d vec2 = neigh_pos2 - m_pos;
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec1));
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec2));

                        angleSum += Vector3d.VectorAngle(normal, vec1);
                        angleSum += Vector3d.VectorAngle(normal, vec2);
                        distSum += vec1.Length;
                        distSum += vec2.Length;

                        srfNodes[i][j].SetNeighPos(vec1);
                        srfNodes[i][j].SetNeighPos(vec2);

                        //Set U_Curvature
                        Plane plane = new Plane(m_pos, normal);
                        if (Math.Abs(plane.DistanceTo(neigh_pos1)) >= 0.001 || Math.Abs(plane.DistanceTo(neigh_pos2)) >= 0.001)  //Check if neighbour points are on same plane, avoid infinite rad.
                        {
                            Circle cir = new Circle(neigh_pos2, m_pos, neigh_pos1);
                            
                            srfNodes[i][j].SetUCurvature(cir.Radius);
                        }
                    }

                    if (j == 0)  //Edge Points @ v = 0
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = srfNodes[i][j + 1].GetPosition();
                        Vector3d vec = neigh_pos - m_pos;
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec));

                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        srfNodes[i][j].SetNeighPos(vec);

                    }
                    else if (j == srfNodes[i].Length - 1) //Edge Points @ v = 1
                    {
                        neighbourCount += 1;
                        Point3d neigh_pos = srfNodes[i][j - 1].GetPosition();
                        Vector3d vec = neigh_pos - m_pos;
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec));

                        angleSum += Vector3d.VectorAngle(normal, vec);
                        distSum += vec.Length;

                        srfNodes[i][j].SetNeighPos(vec);
                    }

                    else  //middle nodes in v direction
                    {
                        neighbourCount += 2;
                        Point3d neigh_pos1 = srfNodes[i][j + 1].GetPosition();
                        Point3d neigh_pos2 = srfNodes[i][j - 1].GetPosition();

                        Vector3d vec1 = neigh_pos1 - m_pos;
                        Vector3d vec2 = neigh_pos2 - m_pos;
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec1));
                        srfNodes[i][j].AddNeighbourAngles(Vector3d.VectorAngle(normal, vec2));

                        angleSum += Vector3d.VectorAngle(normal, vec1);
                        angleSum += Vector3d.VectorAngle(normal, vec2);
                        distSum += vec1.Length;
                        distSum += vec2.Length;

                        srfNodes[i][j].SetNeighPos(vec1);
                        srfNodes[i][j].SetNeighPos(vec2);

                        //Set V_Curvature
                        Plane plane = new Plane(m_pos, normal);
                        if (Math.Abs(plane.DistanceTo(neigh_pos1)) >= 0.001 || Math.Abs(plane.DistanceTo(neigh_pos2)) >= 0.001)  //Check if neighbour points are on same plane, avoid infinite rad.
                        {
                            Circle cir = new Circle(neigh_pos2, m_pos, neigh_pos1);
                            srfNodes[i][j].SetVCurvature(cir.Radius);
                        }
                    }

                    srfNodes[i][j].SetAngleSum(angleSum);
                    srfNodes[i][j].SetNeighbourCount(neighbourCount);
                    srfNodes[i][j].SetDistanceSum(distSum);
                    srfNodes[i][j].Set_UIndex(i);
                    srfNodes[i][j].Set_VIndex(j);
                }
            }
        }

        private void CreateMesh()
        {
            mesh = new Mesh();
            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    Node node = srfNodes[i][j];
                    mesh.Vertices.Add(node.GetPosition());
                }
            }

            for (int i = 0; i < srfNodes.Length - 1; i++)
            {
                for (int j = 0; j < srfNodes[i].Length - 1; j++)
                {
                    int a = (i * srfNodes[i].Length) + j;
                    int b = (i * srfNodes[i].Length) + (j + 1);
                    int c = ((i + 1) * srfNodes[i].Length) + j + 1;
                    int d = ((i + 1) * srfNodes[i].Length) + j;

                    mesh.Faces.AddFace(a, b, c, d);
                }
            }
            mesh.Normals.ComputeNormals();
            mesh.Compact();
        }


        //Public Methods

        /// <summary>
        /// This methods sets the index of nodes as the panel number which holds said nodes.
        /// </summary>
        public void Set_NodePanelIndices(int panel_index)
        {
            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    srfNodes[i][j].Set_PanelIndex(panel_index);
                }
            }
        }
        public Node[][] GetNodes()
        {
            return srfNodes;
        }

        public List<Point3d> GetSurfacePoints()
        {
            List<Point3d> pts = new List<Point3d>();

            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    pts.Add(srfNodes[i][j].GetPosition());
                }
            }
            return pts;
        }
        public List<double> Get_UCurvatures()
        {
            List<double> vals = new List<double>();
            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    vals.Add(srfNodes[i][j].GetUCurvature());
                }
            }
            return vals;
        }
        public List<double> Get_VCurvatures()
        {
            List<double> vals = new List<double>();
            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    vals.Add(srfNodes[i][j].GetVCurvature());
                }
            }
            return vals;
        }
        public Mesh GetSurfaceMesh()
        {
            return mesh;
        }
        public void MoveSurfaceNodes(Vector3d motion)
        {
            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    srfNodes[i][j].MoveNode(motion);
                }
            }

            Transform transform = Transform.Translation(motion);
            this.mesh.Transform(transform);
        }

        public void AdjustPanelRotation()
        {
            Rhino.Geometry.Collections.MeshVertexNormalList normals = mesh.Normals;

            //Vector3d vec = new Vector3d();
            //Transform t = Transform.Identity;
            //if (uDivs > vDivs)
            //{
            //    vec = srfNodes[srfNodes.Length-1][0].m_pos - srfNodes[0][0].m_pos;
            //    t = Transform.Rotation(vec, new Vector3d(1, 0, 0), srfNodes[0][0].m_pos);
            //}
            //else
            //{
            //    vec = srfNodes[0][srfNodes[0].Length-1].m_pos - srfNodes[0][0].m_pos;
            //    t = Transform.Rotation(vec, new Vector3d(0, 1, 0), srfNodes[0][0].m_pos);
            //}
            //Vector3d massAdd = new Vector3d();
            //for (int i = 0; i < normals.Count; i++)
            //{
            //    massAdd += normals[i];
            //}
            //massAdd.Unitize();
            

            Plane plane = new Plane(srfNodes[0][0].m_pos, srf.PointAt(1, 0), srf.PointAt(0, 1));
            Plane plane_To = new Plane(srfNodes[0][0].m_pos, Vector3d.XAxis, Vector3d.YAxis);
            Transform t = Transform.PlaneToPlane(plane, plane_To);

            
            mesh.Transform(t);

            for (int i = 0; i < srfNodes.Length; i++)
            {
                for (int j = 0; j < srfNodes[i].Length; j++)
                {
                    srfNodes[i][j].m_pos.Transform(t);
                }
            }
        }
    }

    ///Class storing Node data wrt neighbours.
    public class Node
    {
        public Point3d m_pos;
        private int m_neighbourCount;
        private double m_angleSum;
        private double m_distanceSum;
        private List<Vector3d> m_neighPos; //All neighbour positions wrt NodePosition
        private List<double> m_neighAngles; //angle wrt every neightbour
        private double uCurvature;
        private double vCurvature;

        private int m_panelIndex;
        private int m_uIndex;
        private int m_vIndex;

        //Node as Particle
        public Vector3d acceleration;
        public  Vector3d velocity;
        private static double drag = 0.5;
        private static double timeStep = 0.05;


        public Node()
        {
            this.m_pos = new Point3d();
            this.m_neighPos = new List<Vector3d>();
            this.m_neighAngles = new List<double>();

            this.acceleration = new Vector3d();
            this.velocity = new Vector3d();
        }

        public Node(Point3d position)
        {
            this.m_pos = position;
            this.m_neighAngles = new List<double>();
            this.m_neighPos = new List<Vector3d>();

            this.velocity = new Vector3d();
            this.acceleration = new Vector3d();
        }

        public void MoveNode(Vector3d motion)
        {
            Transform transform = Transform.Translation(motion);
            this.m_pos.Transform(transform);
        }

        /// <summary>
        /// To use in PinBed, as pins can be adjusted based on vectors. (Particle System)
        /// </summary>
        public void AdjustNode()
        {
            acceleration -= velocity * drag;
            velocity += acceleration * timeStep;
            this.m_pos += velocity * timeStep;
            acceleration *= 0;

        }
        /// <summary>
        /// Add different forces acting on node
        /// </summary>
        public void AddForce(Vector3d force)
        {
            acceleration += force;
        }

        //Set Methods
        public void SetPosition(Point3d pos)
        {
            this.m_pos = pos;
        }
        public void SetNeighbourCount(int count)
        {
            this.m_neighbourCount = count;
        }
        public void SetAngleSum(double angleSum)
        {
            this.m_angleSum = angleSum;
        }
        public void SetDistanceSum(double totalDistances)
        {
            this.m_distanceSum = totalDistances;
        }
        public void SetNeighPos(Vector3d pos)
        {
            m_neighPos.Add(pos);
        }
        public void AddNeighbourAngles(double angle)
        {
            this.m_neighAngles.Add(angle);
        }

        public void SetUCurvature(double value)
        {
            this.uCurvature = value;
        }
        public void SetVCurvature(double value)
        {
            this.vCurvature = value;
        }
        public void Set_PanelIndex(int index)
        {
            this.m_panelIndex = index;
        }
        public void Set_UIndex(int index)
        {
            this.m_uIndex = index;
        }
        public void Set_VIndex(int index)
        {
            this.m_vIndex = index;
        }



        //Get methods
        public Point3d GetPosition()
        {
            return m_pos;
        }
        public int GetNeighbourCount()
        {
            return m_neighbourCount;
        }
        public double GetAngleSum()
        {
            return m_angleSum;
        }
        public double GetDistanceSum()
        {
            return m_distanceSum;
        }
        public List<Vector3d> GetNeighPos()
        {
            return m_neighPos;
        }
        public List<double> GetNeighbourAngles()
        {
            return m_neighAngles;
        }
        public double GetUCurvature()
        {
            return uCurvature;
        }
        public double GetVCurvature()
        {
            return vCurvature;
        }

        public int Get_PanelIndex()
        {
            return this.m_panelIndex;
        }
        public int Get_UIndex()
        {
            return this.m_uIndex;
        }
        public int Get_VIndex()
        {
            return this.m_vIndex;
        }
    }

    public class Panel
    {
        private PinBed bed;
        private Surface srf;
        


        private List<SurfaceNodes> panels;
        private List<Surface> divisions;

        private bool extendedU;
        private bool extendedV;

        private double srf_uLength;
        private double srf_vLength;
        private double panel_uSize;
        private double panel_vSize;

        public int uDivs;
        public int vDivs;

        public Panel(PinBed bed, Surface srf)
        {
            this.srf = srf;
            this.bed = bed;
            this.panels = new List<SurfaceNodes>();
            this.divisions = new List<Surface>();
            

            this.extendedU = false;
            this.extendedV = false;

            Interval uInterval = srf.Domain(0);
            Interval vInterval = srf.Domain(1);
            this.srf_uLength = uInterval.T1-uInterval.T0; //Surface length in U direction
            this.srf_vLength = vInterval.T1-vInterval.T0; //Surface length in V direction

            Panelize();
        }

        private void Panelize()
        {
            this.panel_uSize = bed.uDiv * bed.uSize;
            this.panel_vSize = bed.vDiv * bed.vSize;
            this.uDivs = (int)(srf_uLength / (panel_uSize)); //Number of panels in u direction
            this.vDivs = (int)(srf_vLength / (panel_vSize)); // Number of panels in v direction


            if (uDivs * panel_uSize < srf_uLength)  //If panels sum are smaller than surface size, than increase surface size and panelsize.
            {
                srf.Extend(0, new Interval(0, (uDivs + 1) * panel_uSize));
                this.extendedU = true;
                this.uDivs += 1;  //Create one more layer of panels in u direction.
            }
            if (vDivs * panel_vSize < srf_vLength)
            {
                srf.Extend(1, new Interval(0, (vDivs + 1) * panel_vSize));
                this.extendedV = true;
                this.vDivs += 1; //Create one more layer of panels in v direction.
            }

            for (int i = 0; i <= uDivs; i++)
            {
                for (int j = 0; j <= vDivs; j++)
                {
                    if (i * panel_uSize < srf_uLength && j * panel_vSize < srf_vLength)
                    {
                        Interval uInt = new Interval(i * panel_uSize, (i + 1) * panel_uSize);
                        Interval vInt = new Interval(j * panel_vSize, (j + 1) * panel_vSize);

                        Surface temp = srf.Trim(uInt, vInt);
                        this.divisions.Add(temp);
                    }
                }
            }
            MakeSurfaceNodes();
            Set_NodeIndices();
            //MovePanels();
        }

        private void MakeSurfaceNodes()
        {
            for (int i = 0; i < divisions.Count; i++)
            {
                SurfaceNodes temp = new SurfaceNodes(divisions[i], bed.uDiv, bed.vDiv);
                this.panels.Add(temp);
            }
        }

        private void Set_NodeIndices()
        {
            for (int i = 0; i < panels.Count; i++)
            {
                panels[i].Set_NodePanelIndices(i);
            }
        }
        private void MovePanels()
        {
            for (int i = 0; i < panels.Count; i++)
            {
                SurfaceNodes curr_srfNodes = panels[i];

                List<Point3d> pts = curr_srfNodes.GetSurfacePoints();

                Point3d lowest = new Point3d();
                double minZ = double.MaxValue;
                for (int j = 0; j < pts.Count; j++)
                {
                    if (pts[j].Z < minZ)
                    {
                        minZ = pts[j].Z;
                        lowest = pts[j];
                    }
                }

                Point3d projection = new Point3d(lowest.X, lowest.Y, 0); //get projection on XY plane
                Vector3d motion = projection - lowest; //Move all nodes as per this vector

                int xMotion = (int)(i / vDivs);
                int yMotion = (int)(i % vDivs);

                Vector3d vec = new Vector3d(xMotion * (4 * bed.uSize), yMotion * (4 * bed.vSize), 0);
                vec += motion;

                panels[i].MoveSurfaceNodes(vec);
            }
        }

        public List<SurfaceNodes> GetPanels_SurfaceNodes()
        {
            return panels;
        }
        public List<Mesh> GetPanels_Mesh()
        {
            List<Mesh> meshes = new List<Mesh>();

            for (int i = 0; i < panels.Count; i++)
            {
                meshes.Add(panels[i].GetSurfaceMesh());
            }

            return meshes;
        }
        public List<Surface> GetPanels_Surface()
        {
            return divisions;
        }
        public List<bool> CheckExtension()
        {
            List<bool> bools = new List<bool>();
            bools.Add(this.extendedU);
            bools.Add(this.extendedV);

            return bools;
        } //Returns 2 bools if surface is extended in U or V direction
        public List<int> GetPanelCount()  //Returns list of 2 integers
        {
            List<int> list = new List<int>();

            list.Add(uDivs);
            list.Add(vDivs);

            return list;
        }

        public void Optimize_PanelOrientation()
        {
            for (int i = 0; i < panels.Count; i++)
            {
                panels[i].AdjustPanelRotation();
            }
        }
    }  
}
