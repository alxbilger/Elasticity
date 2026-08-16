#include <Binding_HyperelasticMaterial.h>
#include <Elasticity/component/HyperelasticMaterial.h>
#include <SofaPython3/PythonEnvironment.h>
#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>
#include <string>

namespace elasticity::python3
{
namespace py { using namespace pybind11; }
using namespace std::string_literals;

/**
 * Trampoline class for HyperelasticMaterial
 */
template <class TDataTypes>
class PyHyperelasticMaterial : public HyperelasticMaterial<TDataTypes>
{
public:
    SOFA_CLASS(PyHyperelasticMaterial, HyperelasticMaterial<TDataTypes>);
    using HyperelasticMaterial<TDataTypes>::HyperelasticMaterial;
    using StressTensor = HyperelasticMaterial<TDataTypes>::StressTensor;
    using TangentModulus = HyperelasticMaterial<TDataTypes>::TangentModulus;

    void init() override;
    StressTensor firstPiolaKirchhoffStress(Strain<TDataTypes>& strain) override;
    TangentModulus materialTangentModulus(Strain<TDataTypes>& strain) override;
};

template <class TDataTypes>
void PyHyperelasticMaterial<TDataTypes>::init()
{
    HyperelasticMaterial<TDataTypes>::init();
    sofapython3::PythonEnvironment::gil acquire;
    PYBIND11_OVERRIDE(void, HyperelasticMaterial<TDataTypes>, init, );
}

template <class TDataTypes>
auto PyHyperelasticMaterial<TDataTypes>::firstPiolaKirchhoffStress(Strain<TDataTypes>& strain) -> StressTensor
{
    const auto& F = strain.deformationGradient();
    PYBIND11_OVERLOAD_PURE(StressTensor, HyperelasticMaterial<TDataTypes>, firstPiolaKirchhoffStress, F);
}

template <class TDataTypes>
auto PyHyperelasticMaterial<TDataTypes>::materialTangentModulus(Strain<TDataTypes>& strain) -> TangentModulus
{
    using Real = sofa::Real_t<TDataTypes>;
    using Callable = std::function<Real(sofa::Size, sofa::Size, sofa::Size, sofa::Size)>;

    sofapython3::PythonEnvironment::gil acquire;

    py::function override = py::get_override(
        static_cast<const HyperelasticMaterial<TDataTypes>*>(this),
        "materialTangentModulus");

    if (!override)
    {
        pybind11::pybind11_fail(
            "Tried to call pure virtual function "
            "\"HyperelasticMaterial::materialTangentModulus\"");
    }

    const auto& F = strain.deformationGradient();
    py::object pythonCallable = override(F);

    if (!PyCallable_Check(pythonCallable.ptr()))
    {
        throw py::type_error("materialTangentModulus must return a callable accepting 4 parameters");
    }

    auto callable = [pythonCallable](sofa::Size i, sofa::Size j, sofa::Size k, sofa::Size l) -> Real
    {
        return pythonCallable(i, j, k, l).template cast<Real>();
    };

    return TangentModulus(callable);
}

template <class TDataTypes>
void declareHyperElasticMaterial(pybind11::module &m)
{
    const std::string pyclassName = "HyperElasticMaterial"s + TDataTypes::Name();

    using Class = HyperelasticMaterial<TDataTypes>;
    using Trampoline = PyHyperelasticMaterial<TDataTypes>;
    using Base = sofa::core::objectmodel::BaseComponent;

    py::class_<
        Class,
        Base,
        Trampoline,
        sofapython3::py_shared_ptr<Class>
    > c(m, pyclassName.c_str(), py::dynamic_attr());

    // Python constructor
    c.def(py::init([](py::args &args, py::kwargs &kwargs)
    {
        auto material = sofa::core::sptr<PyHyperelasticMaterial<TDataTypes>> (new PyHyperelasticMaterial<TDataTypes>());
        material->f_listening.setValue(true);

        if (args.size() == 1)
            material->setName(py::cast<std::string>(args[0]));

        py::object cc = py::cast(material);
        for (auto kv : kwargs)
        {
            std::string key = py::cast<std::string>(kv.first);
            py::object value = py::reinterpret_borrow<py::object>(kv.second);
            if (key == "name")
            {
                if (!args.empty())
                {
                    throw py::type_error("The name is set twice as a "
                                        "named argument='" + py::cast<std::string>(value) + "' and as a"
                                                                                            "positional argument='" +
                                        py::cast<std::string>(args[0]) + "'.");
                }
            }
            sofapython3::BindingBase::SetAttr(cc, key, value);
        }
        return material;
    }));

    sofapython3::PythonFactory::registerType<HyperelasticMaterial<TDataTypes>>([](sofa::core::objectmodel::Base* object)
    {
        return py::cast(dynamic_cast<HyperelasticMaterial<TDataTypes>*>(object));
    });
}

void moduleAddHyperelasticMaterial(pybind11::module& m)
{
    declareHyperElasticMaterial<sofa::defaulttype::Vec2Types>(m);
    declareHyperElasticMaterial<sofa::defaulttype::Vec3Types>(m);
}

}  // namespace elasticity::python3
