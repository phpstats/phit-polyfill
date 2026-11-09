<?php

declare(strict_types=1);

use Phit\Polyfill\Internal\Jit;

it('passes the shared IV references and public fit contract', function () {
    if (Jit::compiling()) $this->markTestSkipped('the Phit functions refuse JIT compilation');
    require __DIR__.'/../reference/iv.php';
    expect($nullSw['first_stage']['a']['sw_f_statistic'])->toBeNull();
});
