<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Jit;

it('matches sandwich covariance and inference for weighted OLS and IV', function () {
    if (Jit::compiling()) $this->markTestSkipped('the Phit functions refuse JIT compilation');
    require __DIR__.'/../reference/vcov.php';
    expect($cluster['inference_df'])->toBe(3);
});
