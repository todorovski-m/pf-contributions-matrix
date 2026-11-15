function [ok, msg] = pfc_expect_error(fun, id)
%PFC_EXPECT_ERROR  Call FUN (a function handle) and check that it raises an error with identifier ID.
ok = false;
try
    fun();
    msg = 'no error was raised';
catch err
    ok = strcmp(err.identifier, id);
    msg = sprintf('identifier ''%s'' (expected ''%s''): %s', err.identifier, id, err.message);
end
end
