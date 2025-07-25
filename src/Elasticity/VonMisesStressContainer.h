#pragma once
#include <sofa/config.h>

namespace elasticity
{

template<class Real>
struct VonMisesStressContainer
{
    void resize(sofa::Size numberOfNodes)
    {
        m_counter.resize(numberOfNodes);
        m_stressValues.resize(numberOfNodes);
    }

    void clear()
    {
        m_counter.clear();
        m_stressValues.clear();
    }

    void addVonMisesStress(sofa::Index nodeId, Real vonMisesStress)
    {
        m_stressValues[nodeId] += vonMisesStress;
        ++m_counter[nodeId];
    }

    void accept()
    {
        for (sofa::Size i = 0; i < m_stressValues.size(); ++i)
        {
            if (m_counter[i] != 0)
            {
                m_stressValues[i] /= m_counter[i];
            }
        }
    }

    const sofa::type::vector<Real>& getStressValues() const
    {
        return m_stressValues;
    }

private:
    sofa::type::vector<sofa::Size> m_counter;
    sofa::type::vector<Real> m_stressValues;
};
}
