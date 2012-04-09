
#ifndef Y_MESHTYPES_H
#define Y_MESHTYPES_H

#include <core_api/object3d.h>
#include <yafraycore/triangle.h>


__BEGIN_YAFRAY


struct uv_t
{
	uv_t(GFLOAT _u, GFLOAT _v): u(_u), v(_v) {};
	GFLOAT u, v;
};

class triangle_t;
class vTriangle_t;
class triangleInstance_t;
class triangleObjectInstance_t;

/*!	meshObject_t holds various polygonal primitives
*/

class YAFRAYCORE_EXPORT meshObject_t: public object3d_t
{
	friend class vTriangle_t;
	friend class bsTriangle_t;
	friend class scene_t;
	public:
		meshObject_t(int ntris, bool hasUV=false, bool hasOrco=false);
		/*! the number of primitives the object holds. Primitive is an element
			that by definition can perform ray-triangle intersection */
		int numPrimitives() const { return triangles.size() + s_triangles.size(); }
		int getPrimitives(const primitive_t **prims) const;
		
		primitive_t* addTriangle(const vTriangle_t &t);
		primitive_t* addBsTriangle(const bsTriangle_t &t);
		
		//void setContext(std::vector<point3d_t>::iterator p, std::vector<normal_t>::iterator n);
		void setLight(const light_t *l){ light=l; }
		void finish();
	protected:
		std::vector<vTriangle_t> triangles;
		std::vector<bsTriangle_t> s_triangles;
		std::vector<point3d_t> points;
		std::vector<normal_t> normals;
		std::vector<int> uv_offsets;
		std::vector<uv_t> uv_values;
		bool has_orco;
		bool has_uv;
		bool has_vcol;
		bool is_smooth;
		const light_t *light;
};

/*!	This is a special version of meshObject_t!
	The only difference is that it returns a triangle_t instead of vTriangle_t,
	see declaration if triangle_t for more details!
*/

class YAFRAYCORE_EXPORT triangleObject_t: public object3d_t
{
    friend class triangle_t;
    friend class triangleInstance_t;
	friend class scene_t;
	friend class triangleObjectInstance_t;
	public:
		triangleObject_t() : has_orco(false), has_uv(false), is_smooth(false), normals_exported(false) { /* Empty */ }
		triangleObject_t(int ntris, bool hasUV=false, bool hasOrco=false);
		/*! the number of primitives the object holds. Primitive is an element
			that by definition can perform ray-triangle intersection */
		virtual int numPrimitives() const { return triangles.size(); }
		virtual int getPrimitives(const triangle_t **prims);
		
		triangle_t* addTriangle(const triangle_t &t);
                int add_triangle_with_uv_indices(const triangle_t &t, int const uv0, int const uv1, int const uv2);

                int add_point(point3d_t p) { points.push_back(p); return (points.size() - 1); }
                // int add_normal(point3d_t p) { points.push_back(p); return (points.size() - 1); }

		virtual void finish();

        inline virtual vector3d_t getVertexNormal(int index) const
        {
            return vector3d_t(normals[index]);
        }

        inline virtual point3d_t getVertex(int index) const
        {
            return points[index];
        }

        inline std::vector<triangle_t> const& getTriangles() const
        {
            return triangles;
        }

        inline void setTriangles(std::vector<triangle_t> const& new_triangles)
        {
            triangles = new_triangles;
        }

        void set_use_for_pbgi(bool use_pbgi)
        {
            use_for_pbgi = use_pbgi;
        }

        bool do_use_for_pbgi() const
        {
            return use_for_pbgi;
        }

        triangleObject_t(triangleObject_t const& obj)
        {
            triangles = obj.triangles;
            points = obj.points;
            normals = obj.normals;
            uv_offsets = obj.uv_offsets;
            uv_values = obj.uv_values;
            has_orco = obj.has_orco;
            has_uv = obj.has_uv;
            is_smooth = obj.is_smooth;
            normals_exported = obj.normals_exported;
            use_for_pbgi = obj.use_for_pbgi;

            for (std::size_t i = 0; i < triangles.size(); ++i)
            {
                triangles[i].setObject(this);
            }
        }


	private:
        std::vector<triangle_t> triangles;
		std::vector<point3d_t> points;
		std::vector<normal_t> normals;
		std::vector<int> uv_offsets;
		std::vector<uv_t> uv_values;
	protected:
		bool has_orco;
		bool has_uv;
		bool is_smooth;
		bool normals_exported;
                bool use_for_pbgi;
};

class YAFRAYCORE_EXPORT triangleObjectInstance_t: public triangleObject_t
{
    friend class triangleInstance_t;
	friend class scene_t;
	public:
		triangleObjectInstance_t(triangleObject_t *base, matrix4x4_t obj2World);
		/*! the number of primitives the object holds. Primitive is an element
			that by definition can perform ray-triangle intersection */
		virtual int numPrimitives() const { return triangles.size(); }
		virtual int getPrimitives(const triangle_t **prims);
		
		virtual void finish();

        inline virtual vector3d_t getVertexNormal(int index) const
        {
            return vector3d_t(objToWorld * mBase->normals[index]);
        }

        inline virtual point3d_t getVertex(int index) const
        {
            return objToWorld * mBase->points[index];
        }

	private:
        std::vector<triangleInstance_t> triangles;
        matrix4x4_t objToWorld;
        triangleObject_t* mBase;
};

#include <yafraycore/triangle_inline.h>

__END_YAFRAY


#endif //Y_MESHTYPES_H

