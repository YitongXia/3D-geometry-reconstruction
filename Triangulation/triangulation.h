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

#ifndef TRIANGULATION_H
#define TRIANGULATION_H

#include <easy3d/viewer/viewer.h>

#include "./vector.h"
#include "./matrix.h"

namespace easy3d {
    class Texture;
    class TrianglesDrawable;
    class LinesDrawable;
}

class Triangulation : public easy3d::Viewer {
public:
    Triangulation(
            const std::string &title,
            const std::string &image_point_file_0,
            const std::string &image_point_file_1
    );

protected:
    std::string usage() const override;

    bool key_press_event(int key, int modifiers) override;

    bool triangulation(
            double fx, double fy,     /// input: the focal lengths (same for both cameras)
            double cx, double cy,     /// input: the principal point (same for both cameras)
            const std::vector<easy3d::Vector2D> &image_0_points,    /// input: 2D image points in 1st image.
            const std::vector<easy3d::Vector2D> &image_1_points,    /// input: 2D image points in 2nd image.
            std::vector<easy3d::Vector3D> &points_3d,               /// output: reconstructed 3D points
            easy3d::Matrix33 &R,   /// output: recovered rotation of 2nd camera
            easy3d::Vector3D &t    /// output: recovered translation of 2nd camera
    ) const;

    void post_draw() override;
    void cleanup() override;
    void update_model(const std::vector<easy3d::Vector3D> &points);
    void update_image_plane(const easy3d::mat3 &R, const easy3d::vec3 &t);

private:
    std::vector<easy3d::Vector2D> image_0_points_;
    std::vector<easy3d::Vector2D> image_1_points_;
    easy3d::Texture *texture_0_;
    easy3d::Texture *texture_1_;

    easy3d::TrianglesDrawable *image_plane_;
    easy3d::LinesDrawable *view_frustum_;

    bool show_images_;
};


#endif // TRIANGULATION_H
