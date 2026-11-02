<?php

declare(strict_types=1);

namespace Phit;

/** The fit object of a result whose `model` key is "linear". */
final class LinearFit implements Fit
{
    use \Phit\Polyfill\Objects\FitMethods;
}
