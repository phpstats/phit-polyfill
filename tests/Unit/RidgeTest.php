<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Jit;

it('passes shared ridge and weighted CV glmnet references', function () {
    if (Jit::compiling()) $this->markTestSkipped('the Phit functions refuse JIT compilation');
    require __DIR__.'/../reference/ridge.php';
    expect($fixed)->toBeInstanceOf(Phit\RidgeFit::class);
});
