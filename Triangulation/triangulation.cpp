/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"

#include <easy3d/viewer/drawable_points.h>
#include <easy3d/viewer/camera.h>
#include <easy3d/viewer/texture.h>
#include <easy3d/viewer/primitives.h>
#include <easy3d/viewer/drawable_triangles.h>
#include <easy3d/viewer/drawable_lines.h>
#include <easy3d/fileio/resources.h>
#include <easy3d/fileio/point_cloud_io.h>

#include <3rd_party/glfw/include/GLFW/glfw3.h>    // for the KEYs


using namespace easy3d;


Triangulation::Triangulation(const std::string &title,
                             const std::string &image_point_file_0,
                             const std::string &image_point_file_1)
        : Viewer(title), texture_0_(nullptr), texture_1_(nullptr), image_plane_(nullptr), view_frustum_(nullptr),
          show_images_(true)
{
    auto load_image_points = [](const std::string image_point_file) -> std::vector<Vector2D> {
        PointCloud cloud;
        if (easy3d::load_xyz(image_point_file, cloud)) {
            std::vector<Vector2D> points(cloud.n_vertices());
            const std::vector<vec3> &pts = cloud.points();
            for (std::size_t i = 0; i < pts.size(); ++i) {
                points[i] = Vector(2, pts[i].data());
            }
            return points;
        }
        return std::vector<Vector2D>();
    };

    image_0_points_ = load_image_points(image_point_file_0);
    image_1_points_ = load_image_points(image_point_file_1);

    if (image_0_points_.empty())
        LOG(ERROR) << "failed loading image points from file: " << image_point_file_0;
    else
        LOG(ERROR) << image_0_points_.size() << " image points loaded from file 0";

    if (image_1_points_.empty())
        LOG(ERROR) << "failed loading image points from file: " << image_point_file_1;
    else
        LOG(ERROR) << image_1_points_.size() << " image points loaded from file 1";

    if (image_0_points_.size() != image_1_points_.size())
        LOG(ERROR) << "numbers of image points do not match (" << image_0_points_.size() << " != "
                   << image_1_points_.size() << ")";

    if (image_0_points_.size() < 8)
        LOG(ERROR) << "no enough point pairs (" << image_0_points_.size() << " loaded, but at least " << 8
                   << " needed)";
}


std::string Triangulation::usage() const {
    return ("\n============================== Usage ====================================\n"
            "\tCorresponding images points should have been loaded (unless you see    \n"
            "\t\tsome error messages). If no errors, you can press the 'space' key    \n"
            "\t\tto run your 'triangulation()' method.                                \n"
            "\tYou need to implement 'triangulation()' in 'triangulation_method.cpp'. \n"
            "-------------------------------------------------------------------------\n"
            "\tIf reconstruction was successful, the recovered 3D points will be      \n"
            "\tvisualized for you. You can use your mouse to interact with the viewer:\n"
            "\t\t- Left button: rotate                                                \n"
            "\t\t- Right button: move                                                 \n"
            "\t\t- Wheel: zoom in/out                                                 \n"
            "\t\t- 'i': show/hide images                                               \n"
            "-------------------------------------------------------------------------\n"
            "\tSome other tools provided by the viewer:                               \n"
            "\t\t- 's': take a snapshot of the visualization                          \n"
            "\t\t- 'Ctrl' + 's': save the reconstructed 3D points into an '.xyz' file \n"
            "\t\t   (required in your final submission)                               \n"
            "-------------------------------------------------------------------------\n");
}


bool Triangulation::key_press_event(int key, int modifiers) {
    if (key == GLFW_KEY_SPACE) {
        if (image_0_points_.size() != image_1_points_.size()) {
            LOG(ERROR) << "numbers of points do not match (" << image_0_points_.size() << " != " << image_1_points_.size() << ")";
            return EXIT_FAILURE;
        }

        if (image_0_points_.size() < 8) {
            LOG(ERROR) << "no enough point pairs (" << image_0_points_.size() << " loaded, but at least " << 8 << " needed)";
            return false;
        }

        const double fx = 1000;  /// TODO: modify fx to see how it affects the reconstruction (e.g., 970, 1010, 1030, 1050)
        const double fy = 1000;  /// TODO: modify fy to see how it affects the reconstruction
        const double cx = 320;   /// TODO: modify cx to see how it affects the reconstruction (e.g., 200, 250, 350, 400)
        const double cy = 240;   /// TODO: modify cy to see how it affects the reconstruction
        LOG(ERROR) << "TODO: tune the parameters (above line " << __LINE__ << " in file 'triangulation.cpp') and see how"
                      "\n\terrors in camera intrinsic parameters affect the final reconstruction. Please include your findings"
                      "\n\tand insights in the report";
        std::vector<Vector3D> points_3d;
        Matrix33 R;
        Vector3D t;
        bool success = triangulation(fx, fy, cx, cy, image_0_points_, image_1_points_, points_3d, R, t);
        if (success) {
            if (points_3d.empty()) {
                LOG(ERROR) << "triangulation() has returned 'true', but no 3D points returned";
                return false;
            }

            update_model(points_3d);

            const Box3 &box = current_model()->bounding_box();
            if (box.diagonal() > epsilon<double>()) {
                // resize the viewer to have the same size as the image
                const int scale = dpi_scaling() >= 2.0 ? 3 : 2;
                resize(cx * scale, cy * scale);
                // ensure the points are entirely within the view frustum
                camera()->setSceneBoundingBox(current_model()->bounding_box());

                mat3 RR(R(0, 0), R(0, 1), R(0, 2),
                        R(1, 0), R(1, 1), R(1, 2),
                        R(2, 0), R(2, 1), R(2, 2)
                );
                vec3 tt(t.data());

                // update view using the recovered R and r
                camera()->set_from_calibration(fx, fy, 0.0, cx, cy, RR, tt, false);
                update_image_plane(RR, tt);
                return true;
            } else {
                LOG(ERROR) << "the reconstructed points has a diagonal length (of its bounding box): " << box.diagonal()
                           << ". This value is too small and the reconstruction might not be correct";
                return false;
            }
        } else {
            LOG(ERROR) << "triangulation failed";
            return false;
        }
    } else if (GLFW_KEY_I == key) {
        show_images_ = !show_images_;
        if (image_plane_)
            image_plane_->set_visible(show_images_);
        if (view_frustum_)
            view_frustum_->set_visible(show_images_);
        update();
        return true;
    } else
        return Viewer::key_press_event(key, modifiers);
}


void Triangulation::post_draw() {
    Viewer::post_draw();

    if (!show_images_)
        return;

    if (!texture_0_)
        texture_0_ = Texture::create(resource::directory() + "/data/image_0.png");
    if (!texture_1_)
        texture_1_ = Texture::create(resource::directory() + "/data/image_1.png");

    if (!texture_0_ || !texture_1_)
        return;

    int w = width() * dpi_scaling();
    int h = height() * dpi_scaling();
    double size_ratio = 0.2;

    // ----------------------------------------------------------------------------------------------------------------

    int tex0_w = texture_0_->width();
    int tex0_h = texture_0_->height();
    double image_as0 = tex0_w / static_cast<double>(tex0_h);
    double viewer_as0 = width() / static_cast<double>(height());
    if (image_as0 < viewer_as0) {// thin
        tex0_h = static_cast<int>(height() * size_ratio);
        tex0_w = static_cast<int>(tex0_h * image_as0);
    } else {
        tex0_w = static_cast<int>(width() * size_ratio);
        tex0_h = static_cast<int>(tex0_w / image_as0);
    }
    const int x_min0 = 20 * dpi_scaling();
    const int x_max0 = (20 + tex0_w) * dpi_scaling();
    const Rect quad0(x_min0, x_max0, 20 * dpi_scaling(), (20 + tex0_h) * dpi_scaling());
    opengl::draw_quad_filled(quad0, texture_0_->id(), w, h, -0.999f);

    // ----------------------------------------------------------------------------------------------------------------

    int tex1_w = texture_1_->width();
    int tex1_h = texture_1_->height();
    double image_as1 = tex1_w / static_cast<double>(tex1_h);
    double viewer_as1 = width() / static_cast<double>(height());
    if (image_as1 < viewer_as1) {// thin
        tex1_h = static_cast<int>(height() * size_ratio);
        tex1_w = static_cast<int>(tex1_h * image_as1);
    } else {
        tex1_w = static_cast<int>(width() * size_ratio);
        tex1_h = static_cast<int>(tex1_w / image_as1);
    }
    const int x_min1 = x_max0 + 20 * dpi_scaling();
    const int x_max1 = x_min1 + tex1_w * dpi_scaling();
    if (!current_model()) {
        const Rect quad1(x_min1, x_max1, 20 * dpi_scaling(), (20 + tex1_h) * dpi_scaling());
        opengl::draw_quad_filled(quad1, texture_1_->id(), w, h, -0.999f);
    }
}


void Triangulation::cleanup() {
    if (texture_0_)
        delete texture_0_;

    if (texture_1_)
        delete texture_1_;

    Viewer::cleanup();
}


void Triangulation::update_model(const std::vector<Vector3D> &points) {
    if (points.empty())
        return;

    PointCloud *cloud = dynamic_cast<PointCloud *>(current_model());
    if (cloud)
        cloud->clear();
    else {
        cloud = new PointCloud;
        cloud->set_name("triangulation.xyz");
        add_model(cloud, false);
    }
    std::vector<vec3> pts(points.size());
    for (std::size_t i = 0; i < points.size(); ++i)
        pts[i] = vec3(points[i].data());

    for (std::size_t i = 0; i < pts.size(); ++i)
        cloud->add_vertex(pts[i]);
    std::cout << "reconstructed model has " << cloud->n_vertices() << " points" << std::endl;

    PointsDrawable *drawable = cloud->points_drawable("vertices");
    if (!drawable) {
        drawable = cloud->add_points_drawable("vertices");
        drawable->set_default_color(vec4(0.8f, 0.3f, 0.4f, 1.0f));
        drawable->set_point_size(5.0f);
        drawable->set_impostor_type(PointsDrawable::SPHERE);
    }
    drawable->update_vertex_buffer(pts);
}


// create an image plane
void Triangulation::update_image_plane(const mat3 &R, const easy3d::vec3 &t) {
    if (!image_plane_) {
        image_plane_ = new TrianglesDrawable("image_plane");
        image_plane_->set_visible(show_images_);
        image_plane_->update_texcoord_buffer({vec2(0, 0), vec2(1, 0), vec2(1, 1), vec2(0, 1)});
        image_plane_->update_index_buffer({0, 1, 3, 1, 2, 3});
        image_plane_->set_distinct_back_color(true);
        add_drawable(image_plane_);

        if (!texture_1_)
            texture_1_ = Texture::create(resource::directory() + "/data/image_1.png");

        if (texture_1_) {
            image_plane_->set_texture(texture_1_);
            image_plane_->set_use_texture(true);
        }
    }

    if (!view_frustum_) {
        view_frustum_ = new LinesDrawable("view_frustum");
        view_frustum_->set_visible(show_images_);
        add_drawable(view_frustum_);
    }

    // Determine the distance between the image plane and the camera
#if 1
    // I want to put it just behind the points
    double dist = 0;
    for (auto &v : current_model()->points()) {
        double d = dot(v - camera()->position(), camera()->viewDirection());
        dist = std::max(dist, d);
    }
    dist *= 1.01f; // push slightly behind the points

    const double fov = camera()->fieldOfView();
    const double hw_ratio = static_cast<double>(height()) / width();
    const double halfHeight = dist * tan(fov * 0.5);
    const double halfWidth = halfHeight / hw_ratio;
#else
    const double fov = camera()->fieldOfView();
    const double hw_ratio = static_cast<double>(height()) / width();
    const double halfWidth = camera()->sceneRadius() * 2;  // just a rought estimation, not accurate at all
    const double halfHeight = halfWidth * hw_ratio;
    const double dist = halfHeight / tan(fov * 0.5);
#endif

    // coordinates in camera frame
    const vec3 c(0.0f, 0.0f, 0.0f);  // camera center
    std::vector<vec3> corners = {
            vec3(-halfWidth, -halfHeight, -dist),
            vec3(halfWidth, -halfHeight, -dist),
            vec3(halfWidth, halfHeight, -dist),
            vec3(-halfWidth, halfHeight, -dist)
    };

    Box3 box = current_model()->bounding_box();
    // convert to world coordinate system
    for (auto &p : corners) {
        p = camera()->worldCoordinatesOf(p);
        box.add_point(p);
    }
    box.add_point(camera()->position());
    camera()->setSceneBoundingBox(box);
    camera()->setPivotPoint(current_model()->bounding_box().center());

    image_plane_->update_vertex_buffer(corners);

    corners.push_back(camera()->position());
    view_frustum_->update_vertex_buffer(corners);
    view_frustum_->update_index_buffer({0, 4, 1, 4, 2, 4, 3, 4});
}