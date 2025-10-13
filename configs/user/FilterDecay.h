#include "ActVAction.h"

namespace ActAlgorithm
{
class FilterDecay : public VAction
{
public:
    double fMinLength {};      //!< Min length in X of the cluster set as reference (pads)

public:
    FilterDecay() : VAction("FilterDecay") {}

    void ReadConfiguration(std::shared_ptr<ActRoot::InputBlock> block) override;
    void Run() override;
    void Print() const override;
};
} // namespace ActAlgorithm
