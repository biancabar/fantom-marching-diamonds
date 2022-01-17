#include <fantom/algorithm.hpp>
#include <fantom/datastructures/interfaces/Field.hpp>
#include <fantom/graphics.hpp>
#include <fantom/register.hpp>
#include <math.h>

#include <fantom-plugins/utils/Graphics/HelperFunctions.hpp>
#include <fantom-plugins/utils/Graphics/ObjectRenderer.hpp>

#include <stdexcept>
#include <vector>

using namespace fantom;

namespace
{

    class MarchingDiamondsAlgorithm : public VisAlgorithm
    {

    public:
        struct Options : public VisAlgorithm::Options
        {
            Options(fantom::Options::Control &control)
                : VisAlgorithm::Options(control)
            {
                add<Field<2, Scalar>>("Field", "A 2D scalar field", definedOn<Grid<2>>(Grid<2>::Points));
                add<Color>("Color", "The color of the graphics.", Color(0.75, 0.75, 0.0));
                add<double>("Threshold", "The threshold which points to show.", 0.0008);
                add<double>("Radius", "The size of the scalar-ellipsoids.", 0.1);
            }
        };

        struct VisOutputs : public VisAlgorithm::VisOutputs
        {
            VisOutputs(fantom::VisOutputs::Control &control)
                : VisAlgorithm::VisOutputs(control)
            {
                addGraphics("Ellipsoids");
            }
        };

        MarchingDiamondsAlgorithm(InitData &data)
            : VisAlgorithm(data)
        {
        }

        virtual void execute(const Algorithm::Options &options, const volatile bool &abortFlag) override
        {
            std::shared_ptr<const Field<2, Scalar>> field = options.get<Field<2, Scalar>>("Field");
            std::shared_ptr<const Function<Scalar>> function = options.get<Function<Scalar>>("Field");
            Color color = options.get<Color>("Color");
            double threshold = options.get<double>("Threshold");
            double radius = options.get<double>("Radius");

            //check if input field is set
            if (!field)
            {
                debugLog() << "Input field not set." << std::endl;
                return;
            }

            std::shared_ptr<const Grid<2>> grid = std::dynamic_pointer_cast<const Grid<2>>(function->domain());

            //check if grid can be used or is correct
            if (!grid)
            {
                throw std::logic_error("Wrong type of grid!");
            }

            //get control points of the grid
            const ValueArray<Point2> &points = grid->points();
            std::vector<Point3> vertices;

            auto evaluator = field->makeEvaluator();

            //iterate through every point of the field, checking the value against the threshold
            for (size_t i = 0; i < grid->numPoints(); i++)
            {
                Point2 point = points[i];
                Point3 pointFlat = {point[0], point[1], (double)0};
                evaluator->reset(point);

                auto value = evaluator->value();

                if (value[0] > threshold)
                {
                    vertices.insert(vertices.end(), pointFlat);
                }
            }

            if (abortFlag)
            {
                return;
            }

            auto const &system = graphics::GraphicsSystem::instance();
            std::string resourcePath = PluginRegistrationService::getInstance().getResourcePath("utils/Graphics");

            auto performanceObjectRenderer = std::make_shared<graphics::ObjectRenderer>(system);
            performanceObjectRenderer->reserve(graphics::ObjectRenderer::ObjectType::SPHERE, vertices.size());

            //add sphere to renderer
            for (size_t i = 0; i < vertices.size(); i++)
            {
                performanceObjectRenderer->addSphere(vertices[i], radius, color);
            }

            setGraphics("Ellipsoids", performanceObjectRenderer->commit());

            /************************************/
            std::vector<Point2> intersect = calcIntersection(1, 0.5,1,0.5,1);
            
            infoLog() << "Intersections: " << intersect[0] << std::endl;
        }

    private:

        /*calc intersections of edge e with bilinear plane (hyperbola for fixed iso-value c) in mapped unit square 
        
                    s3                   
                   /  \                 s1 -> (0,1)         s1----s2
                  /    \                s2 -> (1,1)         |    / |
                s1     s0      --->     s0 -> (1,0)   --->  |   e  |
                  \    /                s3 -> (0,0)         | /    |
                   \  /                                     s3----s0
                    s2                  
        */
        std::vector<Point2> calcIntersection(double c, double s0, double s1, double s2, double s3)
        {
            std::vector<Point2> intersections;
            double denominator = 2*( s0 + s1 - s2 + s3 );
            double result;

            //check if denominator is too close/equals zero, return empty as error
            if (abs(denominator) < 0.001) {
                return intersections;
            }

            double w = 4 * (s3 - c)*(s0 + s1 - s2 + s3) + ((s0 + s1) * (s0 + s1));
            //check if root is real, return empty as error (only real values)
            if (w < 0) {
                return intersections;
            }

            w = sqrt(w);
            if (w == 0) {
                result = (s0 + s1) / denominator;
                intersections.push_back({Point2({result, result})});
            } else {
                result = (s0 + s1 + w) / denominator;
                intersections.push_back({Point2({result, result})});

                result = (s0 + s1 - w) / denominator;
                intersections.push_back({Point2({result, result})});
            }

            return intersections;
        }

        //maps coordinats inside the unit square to the reference diamond R via bilinear transformation
        //ONLY for specific intersection condition where x = y
        Point2 bilinearTransformUTR(double y)
        {
            Point2 v1 = {0, 2};
            Point2 v2 = {0, -2};

            return v1 + (v2 * y * y);
        }
        
    };

    AlgorithmRegister<MarchingDiamondsAlgorithm> dummy("Iso-Surface/MarchingDiamonds2D",
                                                           "Show scalar values in 2D grid over certain threshold.");
} // namespace
